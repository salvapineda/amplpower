import logging
from pathlib import Path

import numpy as np
import pandas as pd
from amplpy import AMPL
from amplpy import add_to_path
from matpowercaseframes import CaseFrames

add_to_path(r"/opt/ampl/")


def compute(args):
    return max(args, key=len)


def array2dict(array):
    """Convert a 2D numpy array to a dictionary."""
    return {(i, j): array[i, j] for i in range(array.shape[0]) for j in range(array.shape[1])}


class PowerSystem:
    def __init__(self, case_file: str):
        print(f"=======Initializing the power system with case file: {case_file}")
        self.case_file = case_file
        self.baseMVA = None
        self.buses = pd.DataFrame()
        self.generators = pd.DataFrame()
        self.branches = pd.DataFrame()
        self.gencost = pd.DataFrame()
        self.nbus = 0
        self.nlin = 0
        self.ngen = 0
        self.epsilon = 1e-6  # Tolerance for zero values

        # Admittance matrices (initialize with zeros)
        self.yff = np.zeros((0, 0))
        self.ytf = np.zeros((0, 0))
        self.yft = np.zeros((0, 0))
        self.ytt = np.zeros((0, 0))

        # Connection matrices
        self.cf = None
        self.ct = None
        self.cg = None  # Matrix for generator connections

        self._load_data()
        self._initialize_matrices()
        self.compute_admittance_matrices()
        self.initialize()
        self.summary()
        self.compute_initial_bigm_dc()
        self.compute_initial_bigm_ac()

    def _load_data(self):
        """Load MATPOWER case data into DataFrames and convert to per unit."""
        try:
            case = CaseFrames(self.case_file)
            # Load data for each component
            self.baseMVA = case.baseMVA
            self.buses = case.bus
            self.buses.reset_index(drop=True, inplace=True)
            self.buses["BUS_I"] -= 1
            self.generators = case.gen
            self.generators.reset_index(drop=True, inplace=True)
            self.generators["GEN_BUS"] -= 1
            self.branches = case.branch
            self.branches.reset_index(drop=True, inplace=True)
            self.branches["F_BUS"] -= 1
            self.branches["T_BUS"] -= 1
            self.gencost = case.gencost
            self.gencost.reset_index(drop=True, inplace=True)
            self.nbus = len(self.buses)  # Number of buses
            self.nlin = len(self.branches)  # Number of branches
            self.ngen = len(self.generators)  # Number of generators

            # Minimum and maximum voltage limits
            self.max_voltage = self.buses["VMAX"].max()
            self.min_voltage = self.buses["VMIN"].min()
            self.max_angle = np.pi / 2
            self.min_angle = -np.pi / 2

            # Convert to per unit
            self.buses["PD"] /= self.baseMVA
            self.buses["QD"] /= self.baseMVA
            self.buses["GS"] /= self.baseMVA
            self.buses["BS"] /= self.baseMVA
            self.generators["PG"] /= self.baseMVA
            self.generators["QG"] /= self.baseMVA
            self.generators["PMAX"] /= self.baseMVA
            self.generators["PMIN"] /= self.baseMVA
            self.generators["QMAX"] /= self.baseMVA
            self.generators["QMIN"] /= self.baseMVA
            self.branches["RATE_A"] /= self.baseMVA
            self.branches["RATE_B"] /= self.baseMVA
            self.branches["RATE_C"] /= self.baseMVA

            # Set default branch limit
            self.default_branch_limit = np.sqrt(self.buses["PD"].sum() ** 2 + self.buses["QD"].sum() ** 2)
            for line_index in range(self.nlin):
                if self.branches.loc[line_index, "RATE_A"] == 0:
                    self.branches.loc[line_index, "RATE_A"] = self.default_branch_limit
        except Exception as e:
            logging.error(f"Error loading data from {self.case_file}: {e}")
            raise

    def _initialize_matrices(self):
        """Initialize matrices for admittance calculations."""
        self.yff = np.zeros(self.nlin, dtype=complex)
        self.ytf = np.zeros(self.nlin, dtype=complex)
        self.yft = np.zeros(self.nlin, dtype=complex)
        self.ytt = np.zeros(self.nlin, dtype=complex)
        self.cf = np.zeros((self.nlin, self.nbus))  # Connection for F_BUS
        self.ct = np.zeros((self.nlin, self.nbus))  # Connection for T_BUS
        self.cg = np.zeros((self.ngen, self.nbus))  # Connection for generators

        # Update generator connection matrix
        for g in range(self.ngen):
            bus = int(self.generators.iloc[g]["GEN_BUS"])  # Ensure index is an integer
            self.cg[g, bus] = 1

    def compute_admittance_matrices(self):
        """Calculate the admittance matrices (yff, ytf, yft, ytt) for the network."""
        for line_index in range(self.nlin):
            branch = self.branches.iloc[line_index]  # Access branch data
            r = branch["BR_R"]
            x = branch["BR_X"]
            b = branch["BR_B"]
            tau = branch["TAP"] if branch["TAP"] != 0 else 1  # Handle TAP=0 case
            theta = branch["SHIFT"]

            # Calculate Y series and shunt admittance
            ys = 1 / (r + 1j * x)

            # Store the admittance components
            self.yff[line_index] = (ys + 1j * 0.5 * b) / (tau**2)
            self.yft[line_index] = -ys / (tau * np.exp(-1j * theta))
            self.ytf[line_index] = -ys / (tau * np.exp(1j * theta))
            self.ytt[line_index] = ys + 1j * 0.5 * b

            # Update bus connection matrices
            f_bus, t_bus = int(branch["F_BUS"]), int(branch["T_BUS"])  # Ensure indices are integers
            self.cf[line_index, f_bus] = 1
            self.ct[line_index, t_bus] = 1

        # Compute additional matrices
        self.yf = np.dot(np.diag(self.yff), self.cf) + np.dot(np.diag(self.yft), self.ct)
        self.yt = np.dot(np.diag(self.ytf), self.cf) + np.dot(np.diag(self.ytt), self.ct)
        self.ysh = self.buses["GS"].values + 1j * self.buses["BS"].values
        self.yb = np.dot(np.transpose(self.cf), self.yf) + np.dot(np.transpose(self.ct), self.yt) + np.diag(self.ysh)

        # Include admittance values in the branch DataFrame
        self.branches["GFF"] = np.real(self.yff)
        self.branches["BFF"] = np.imag(self.yff)
        self.branches["GFT"] = np.real(self.yft)
        self.branches["BFT"] = np.imag(self.yft)
        self.branches["GTF"] = np.real(self.ytf)
        self.branches["BTF"] = np.imag(self.ytf)
        self.branches["GTT"] = np.real(self.ytt)
        self.branches["BTT"] = np.imag(self.ytt)

    def initialize(self, voltages=None, angles=None):
        """Initialize the voltage magnitudes, angles, flows, and generation levels."""
        if voltages is None:
            voltages = np.ones(self.nbus)
        if angles is None:
            angles = np.zeros(self.nbus)
        self.buses["VOL0"] = voltages
        self.buses["ANG0"] = angles
        self.buses["VOLR0"] = voltages * np.cos(angles)
        self.buses["VOLI0"] = voltages * np.sin(angles)

        # Compute flows
        v = voltages * np.exp(1j * angles)
        sf = (self.cf @ v) * np.conj(self.yf @ v)
        st = (self.ct @ v) * np.conj(self.yt @ v)
        self.branches["PF0"] = np.real(sf)
        self.branches["QF0"] = np.imag(sf)
        self.branches["PT0"] = np.real(st)
        self.branches["QT0"] = np.imag(st)

        # Compute generator outputs
        sd = self.buses["PD"].values + 1j * self.buses["QD"].values
        sb = v * np.conj(self.yb @ v)
        sg = sb + sd
        self.generators["PG0"] = np.dot(np.real(sg), self.cg.T)
        self.generators["QG0"] = np.dot(np.imag(sg), self.cg.T)

    def summary(self):
        """Print summary of the network."""
        print(f"Number of buses: {self.nbus}")
        print(f"Number of lines: {self.nlin}")
        print(f"Number of generators: {self.ngen}")
        print(f"baseMVA: {self.baseMVA}")
        print("\nBuses:")
        print(self.buses.head())
        print("\nGenerators:")
        print(self.generators.head())
        print("\nBranches:")
        print(self.branches.head())
        print("\nGenerator Costs:")
        print(self.gencost.head())

    def compute_initial_bigm_dc(self):
        print("=======Computing initial bigM values for DC power flow")
        """Compute Big-M values for the different lines and return them in a DataFrame."""
        self.branches["PFUPDC"] = (1 / self.branches["BR_X"]) * np.pi
        self.branches["PFLODC"] = -(1 / self.branches["BR_X"]) * np.pi

    def compute_initial_bigm_ac(self):
        print("=======Computing initial bigM values for AC power flow")
        """Compute Big-M values for the different lines and return them in a DataFrame."""
        v2max = self.max_voltage**2
        v2min = self.min_voltage**2
        cosmax = self.max_voltage**2
        cosmin = 0
        sinmax = self.max_voltage**2
        sinmin = -(self.max_voltage**2)
        self.branches["PFUPAC"] = self.branches["GFF"] * v2max + self.branches["GFT"] * cosmin + self.branches["BFT"] * sinmax
        self.branches["PFLOAC"] = self.branches["GFF"] * v2min + self.branches["GFT"] * cosmax + self.branches["BFT"] * sinmin
        self.branches["PTUPAC"] = self.branches["GTT"] * v2max + self.branches["GTF"] * cosmin + self.branches["BTF"] * sinmax
        self.branches["PTLOAC"] = self.branches["GTT"] * v2min + self.branches["GTF"] * cosmax + self.branches["BTF"] * sinmin
        self.branches["QFUPAC"] = -self.branches["BFF"] * v2max - self.branches["BFT"] * cosmin + self.branches["GFT"] * sinmin
        self.branches["QFLOAC"] = -self.branches["BFF"] * v2min - self.branches["BFT"] * cosmax + self.branches["GFT"] * sinmax
        self.branches["QTUPAC"] = -self.branches["BTT"] * v2max - self.branches["BTF"] * cosmin + self.branches["GTF"] * sinmin
        self.branches["QTLOAC"] = -self.branches["BTT"] * v2min - self.branches["BTF"] * cosmax + self.branches["GTF"] * sinmax

    def solve_opf(self, opf_type="dc", switching="off", solver="gurobi", time_limit=3600):
        """Solve the optimal power flow problem using AMPL.
        Parameters:
        opf_type (str): Type of optimal power flow ('dc', 'acrect', 'acjabr')
        switching (str): Switching strategy ('off', 'nl', 'bigm')
        solver (str): Solver to use ('gurobi', 'cplex', 'cbc')
        time_limit (int): Time limit for the optimization
        Returns:
        dict: Results of the optimal power flow problem
        """
        # set the status of the lines
        if isinstance(switching, np.ndarray):
            self.branches["BR_STATUS"] = switching
        elif switching == "off":
            self.branches["BR_STATUS"] = 1
        elif switching == "nl":
            self.branches["BR_STATUS"] = 2
        elif switching == "bigm":
            self.branches["BR_STATUS"] = 3

        print(f"=======Solving OPF with AMPL ({opf_type}) with solver {solver}, time limite {time_limit}")
        ampl = AMPL()
        ampl.reset()
        ampl.read(Path(__file__).parent / "opf.mod")

        ampl.set_data(self.buses, "N")
        ampl.set_data(self.generators, "G")
        ampl.set_data(self.branches, "L")
        ampl.set_data(self.gencost)
        ampl.param["CF"] = array2dict(self.cf)
        ampl.param["CT"] = array2dict(self.ct)
        ampl.param["CG"] = array2dict(self.cg)
        ampl.param["OPF_TYPE"] = opf_type
        ampl.param["BASEMVA"] = self.baseMVA
        ampl.param["MAXVOL"] = self.max_voltage
        ampl.param["MINVOL"] = self.min_voltage
        ampl.param["MAXANGLE"] = self.max_angle
        ampl.param["MINANGLE"] = self.min_angle

        ampl.option["mp_options"] = f"mipgap=0.001 threads=1 outlev=1 timelimit={time_limit}"
        ampl.solve(solver=solver)
        solver_status = ampl.solve_result

        if solver_status == "solved" or solver_status == "limit":
            # Extract solution
            genp_values = ampl.get_variable("genp").get_values().to_pandas().values.flatten()
            genq_values = ampl.get_variable("genq").get_values().to_pandas().values.flatten()
            gen_df = pd.DataFrame({"Pg": genp_values, "Qg": genq_values}, index=ampl.get_variable("genp").get_values().to_pandas().index)

            vol_values = ampl.get_variable("vol").get_values().to_pandas().values.flatten()
            ang_values = ampl.get_variable("ang").get_values().to_pandas().values.flatten()
            bus_df = pd.DataFrame({"Vm": vol_values, "Va": ang_values}, index=ampl.get_variable("vol").get_values().to_pandas().index)

            status_values = ampl.get_variable("status").get_values().to_pandas().values.flatten()
            flowpf_values = ampl.get_variable("flowpf").get_values().to_pandas().values.flatten()
            flowpt_values = ampl.get_variable("flowpt").get_values().to_pandas().values.flatten()
            flowqf_values = ampl.get_variable("flowqf").get_values().to_pandas().values.flatten()
            flowqt_values = ampl.get_variable("flowqt").get_values().to_pandas().values.flatten()
            line_df = pd.DataFrame(
                {"switching": status_values, "Pf": flowpf_values, "Pt": flowpt_values, "Qf": flowqf_values, "Qt": flowqt_values},
                index=ampl.get_variable("status").get_values().to_pandas().index,
            )

            return {
                "obj": ampl.get_objective("total_cost").value(),
                "time": ampl.get_value("_solve_time"),
                "gen": gen_df,
                "bus": bus_df,
                "lin": line_df,
            }

        elif solver_status == "infeasible":
            return None
