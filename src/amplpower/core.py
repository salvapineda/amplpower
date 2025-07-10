import logging
from pathlib import Path

import numpy as np
import pandas as pd
from amplpy import AMPL
from matpowercaseframes import CaseFrames

from amplpower.utils import find_min_max


def compute(args):
    return max(args, key=len)


# TODO: Remove compute function at some point


def array2dict(array):
    """Convert a 2D numpy array to a dictionary."""
    return {(i, j): array[i, j] for i in range(array.shape[0]) for j in range(array.shape[1])}


class PowerSystem:
    """PowerSystem class for solving optimal power flow problems."""

    def __init__(self, case_file: str, max_voltage_angle: float = np.pi / 2, min_voltage_angle: float = -np.pi / 2):
        """Initialize the power system with a MATPOWER case file."""
        print(f"=======Initializing the power system with case file: {case_file}")
        self.case_file = case_file
        self.max_angle = max_voltage_angle
        self.min_angle = min_voltage_angle
        self.load_data()
        self.summary()
        self.compute_matrices()
        self.initialize()
        self.compute_voltage_bounds()
        self.compute_bigm_dc()
        self.compute_bigm_ac()
        self.compute_ptdf()
        self.compute_lodf()

    def load_data(self):
        """Load MATPOWER case data into DataFrames and convert to per unit."""
        try:
            case = CaseFrames(self.case_file)
            # Load data for each component
            self.baseMVA = case.baseMVA
            self.buses = case.bus
            self.buses.reset_index(drop=True, inplace=True)

            # Create a mapping for bus numbers
            self.bus_mapping = {bus: idx for idx, bus in enumerate(self.buses["BUS_I"])}
            self.reverse_bus_mapping = {idx: bus for bus, idx in self.bus_mapping.items()}
            self.buses["BUS_I"] = self.buses["BUS_I"].map(self.bus_mapping)

            self.generators = case.gen
            self.generators.reset_index(drop=True, inplace=True)
            self.generators["GEN_BUS"] = self.generators["GEN_BUS"].map(self.bus_mapping)

            self.branches = case.branch
            self.branches.reset_index(drop=True, inplace=True)
            self.branches["F_BUS"] = self.branches["F_BUS"].map(self.bus_mapping)
            self.branches["T_BUS"] = self.branches["T_BUS"].map(self.bus_mapping)

            self.gencost = case.gencost
            self.gencost.reset_index(drop=True, inplace=True)
            self.nbus = len(self.buses)
            self.nlin = len(self.branches)
            self.ngen = len(self.generators)

            # Add default values for generator costs if not provided
            if "COST_2" not in self.gencost.columns:
                self.gencost["COST_2"] = 0

            # Minimum and maximum limits for voltage angle and real/imaginary voltage variables
            self.buses["AMAX"] = self.max_angle
            self.buses["AMIN"] = self.min_angle

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

            # Set default branch limit if not provided
            self.default_branch_limit = np.sqrt(self.buses["PD"].sum() ** 2 + self.buses["QD"].sum() ** 2)
            for line_index in range(self.nlin):
                if self.branches.loc[line_index, "RATE_A"] == 0:
                    self.branches.loc[line_index, "RATE_A"] = self.default_branch_limit

            # Define PFMAX and PFMIN for branches
            self.branches["PFMAX"] = self.branches["RATE_A"]
            self.branches["PFMIN"] = -self.branches["RATE_A"]
            self.branches["QFMAX"] = self.branches["RATE_A"]
            self.branches["QFMIN"] = -self.branches["RATE_A"]

        except Exception as e:
            logging.error(f"Error loading data from {self.case_file}: {e}")
            raise

    def compute_matrices(self):
        """Calculate the admittance matrices (yff, ytf, yft, ytt) for the network."""
        # Initizalize matrices
        self.yff = np.zeros(self.nlin, dtype=complex)
        self.ytf = np.zeros(self.nlin, dtype=complex)
        self.yft = np.zeros(self.nlin, dtype=complex)
        self.ytt = np.zeros(self.nlin, dtype=complex)
        self.cf = np.zeros((self.nlin, self.nbus))  # Connection for F_BUS
        self.ct = np.zeros((self.nlin, self.nbus))  # Connection for T_BUS
        self.cg = np.zeros((self.ngen, self.nbus))  # Connection for generators
        # Compute admittance matrices
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

        # Compute generator connection matrix
        for g in range(self.ngen):
            bus = int(self.generators.iloc[g]["GEN_BUS"])  # Ensure index is an integer
            self.cg[g, bus] = 1

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
        pg_split = np.zeros(self.ngen)
        qg_split = np.zeros(self.ngen)
        for bus in range(self.nbus):
            gen_indices = self.generators[self.generators["GEN_BUS"] == bus].index
            if len(gen_indices) > 0:
                pmax_total = self.generators.loc[gen_indices, "PMAX"].sum()
                if pmax_total > 0:
                    pg_split[gen_indices] = np.real(sg[bus]) * self.generators.loc[gen_indices, "PMAX"] / pmax_total
                    qg_split[gen_indices] = np.imag(sg[bus]) * self.generators.loc[gen_indices, "PMAX"] / pmax_total
        self.generators["PG0"] = pg_split
        self.generators["QG0"] = qg_split

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

    def compute_voltage_bounds(self):
        """Compute bounds for the real and imaginary parts of voltage."""
        print("=======Computing voltage bounds")
        for bus in range(self.nbus):
            vmin = self.buses.loc[bus, "VMIN"]
            vmax = self.buses.loc[bus, "VMAX"]
            amin = self.buses.loc[bus, "AMIN"]
            amax = self.buses.loc[bus, "AMAX"]
            v_critical = [vmin, vmax]
            a_critical = [amin, amax]
            for ac in [0, np.pi / 2, np.pi, -np.pi / 2, -1 * np.pi]:
                if amin <= ac <= amax:
                    a_critical.append(ac)
            vrmin = float("inf")
            vrmax = float("-inf")
            vimin = float("inf")
            vimax = float("-inf")
            for v in v_critical:
                for a in a_critical:
                    vrmin = min(vrmin, v * np.cos(a))
                    vrmax = max(vrmax, v * np.cos(a))
                    vimin = min(vimin, v * np.sin(a))
                    vimax = max(vimax, v * np.sin(a))
            self.buses.loc[bus, "VRMIN"] = vrmin
            self.buses.loc[bus, "VRMAX"] = vrmax
            self.buses.loc[bus, "VIMIN"] = vimin
            self.buses.loc[bus, "VIMAX"] = vimax

    def compute_bigm_dc(self):
        """Compute Big-M values for DC power flow."""
        print("=======Computing Big-M values for DC power flow")
        weights = self.branches["PFMAX"] * self.branches["BR_X"]
        bound_lp = weights.nlargest(self.nlin - 1).sum()

        self.branches["PFUPDC"] = bound_lp / self.branches["BR_X"]
        self.branches["PFLODC"] = -bound_lp / self.branches["BR_X"]

        self.branches["PFUPDC"] = np.ceil(self.branches["PFUPDC"] * 100) / 100
        self.branches["PFLODC"] = np.floor(self.branches["PFLODC"] * 100) / 100

    def compute_bigm_ac(self):
        """Compute Big-M values for active and reactive power flows."""
        print("=======Computing Big-M values for AC power flow using find_min_max")
        for lin_index in range(self.nlin):
            branch = self.branches.iloc[lin_index]
            f_bus = int(branch["F_BUS"])
            t_bus = int(branch["T_BUS"])
            vmaxf, vminf = self.buses.loc[f_bus, "VMAX"], self.buses.loc[f_bus, "VMIN"]
            vmaxt, vmint = self.buses.loc[t_bus, "VMAX"], self.buses.loc[t_bus, "VMIN"]
            amaxf, aminf = self.buses.loc[f_bus, "AMAX"], self.buses.loc[f_bus, "AMIN"]
            amaxt, amint = self.buses.loc[t_bus, "AMAX"], self.buses.loc[t_bus, "AMIN"]

            # Active power flow at "from" bus
            pf_min, pf_max = find_min_max(
                branch["GFF"], branch["GFT"], branch["BFT"], vminf, vmaxf, vmint, vmaxt, aminf - amaxt, amaxf - amint
            )
            self.branches.loc[lin_index, "PFLOAC"] = np.floor(pf_min * 100) / 100
            self.branches.loc[lin_index, "PFUPAC"] = np.ceil(pf_max * 100) / 100

            # Active power flow at "to" bus
            pt_min, pt_max = find_min_max(
                branch["GTT"], branch["GTF"], branch["BTF"], vmint, vmaxt, vminf, vmaxf, amint - amaxf, amaxt - aminf
            )
            self.branches.loc[lin_index, "PTLOAC"] = np.floor(pt_min * 100) / 100
            self.branches.loc[lin_index, "PTUPAC"] = np.ceil(pt_max * 100) / 100

            # Reactive power flow at "from" bus
            qf_min, qf_max = find_min_max(
                -branch["BFF"], -branch["BFT"], branch["GFT"], vminf, vmaxf, vmint, vmaxt, aminf - amaxt, amaxf - amint
            )
            self.branches.loc[lin_index, "QFLOAC"] = np.floor(qf_min * 100) / 100
            self.branches.loc[lin_index, "QFUPAC"] = np.ceil(qf_max * 100) / 100

            # Reactive power flow at "to" bus
            qt_min, qt_max = find_min_max(
                -branch["BTT"], -branch["BTF"], branch["GTF"], vmint, vmaxt, vminf, vmaxf, amint - amaxf, amaxt - aminf
            )
            self.branches.loc[lin_index, "QTLOAC"] = np.floor(qt_min * 100) / 100
            self.branches.loc[lin_index, "QTUPAC"] = np.ceil(qt_max * 100) / 100

            # Cosine of angle difference
            cos_min, cos_max = find_min_max(0, 1, 0, vminf, vmaxf, vmint, vmaxt, aminf - amaxt, amaxf - amint)
            self.branches.loc[lin_index, "COSFTMAX"] = np.ceil(cos_max * 100) / 100
            self.branches.loc[lin_index, "COSFTMIN"] = np.floor(cos_min * 100) / 100

            # Sine of angle difference
            sin_min, sin_max = find_min_max(0, 0, 1, vminf, vmaxf, vmint, vmaxt, aminf - amaxt, amaxf - amint)
            self.branches.loc[lin_index, "SINFTMAX"] = np.ceil(sin_max * 100) / 100
            self.branches.loc[lin_index, "SINFTMIN"] = np.floor(sin_min * 100) / 100

    def compute_ptdf(self):
        """
        Compute the PTDF (Power Transfer Distribution Factor) matrix and store it in self.ptdf.
        The slack bus is the one with BUS_TYPE == 3 in self.buses dataframe.
        """
        nlin = self.nlin
        nbus = self.nbus
        # Find reference bus (BUS_TYPE == 3)
        ref_bus_idx = self.buses.index[self.buses["BUS_TYPE"] == 3][0]
        ref_bus_pos = list(self.buses.index).index(ref_bus_idx)
        # Build branch-to-node incidence matrix A (nlin x nbus)
        A = self.cf - self.ct  # shape (nlin, nbus)
        # Remove slack bus column
        keep = [i for i in range(nbus) if i != ref_bus_pos]
        A_red = A[:, keep]  # shape (nlin, nbus-1)
        # Diagonal matrix of line reactances
        X = np.diag(self.branches["BR_X"].values)  # shape (nlin, nlin)
        Xinv = np.linalg.inv(X)
        # Compute B = A_red.T @ Xinv @ A_red
        B = A_red.T @ Xinv @ A_red  # shape (nbus-1, nbus-1)
        B_inv = np.linalg.inv(B)
        # Compute PTDF: Xinv @ A_red @ B_inv
        PTDF_red = Xinv @ A_red @ B_inv  # shape (nlin, nbus-1)
        # Insert zeros for slack bus column
        PTDF = np.zeros((nlin, nbus))
        PTDF[:, keep] = PTDF_red
        self.ptdf = PTDF

    def compute_lodf(self):
        """
        Compute the LODF (Line Outage Distribution Factor) matrix and store it in self.lodf.
        Uses the PTDF matrix already stored in self.ptdf.
        """
        nlin = self.nlin
        PTDF = self.ptdf  # shape (nlin, nbus)
        # Build branch-to-node incidence matrix A (nlin x nbus)
        A = self.cf - self.ct  # shape (nlin, nbus)
        # H = PTDF @ A.T (A.T is nbus x nlin)
        H = PTDF @ A.T  # shape (nlin, nlin)
        h = np.diag(H)
        denom_mat = np.ones((nlin, nlin)) - np.outer(np.ones(nlin), h)
        LODF = H / denom_mat
        np.fill_diagonal(LODF, 0)
        LODF = LODF - np.eye(nlin)
        self.lodf = LODF

    def set_switching(self, switching):
        """Set the switching status of the branches.

        Explanations:
        - 'off': All lines are connected (no switching, all BR_SWITCH=1).
        - 'nl': Use binary variables for line connection with a non-linear formulation (BR_SWITCH=2).
        - 'bigm': Use binary variables for line connection with a Big-M formulation (BR_SWITCH=3).
        - 'df': Use the values already stored in BR_SWITCH (do not modify).
        - numpy.ndarray: Directly set BR_SWITCH to the provided array.

        Switching statuses:
        0: The line is off.
        1: The line is on.
        2: The line is switchable and modeled with a non-linear approach.
        3: The line is switchable and modeled with a Big-M approach.

        Parameters:
        switching (str or np.ndarray): The switching strategy or array of statuses.
        """
        if isinstance(switching, np.ndarray):
            self.branches["BR_SWITCH"] = switching
        elif switching == "off":
            self.branches["BR_SWITCH"] = 1
        elif switching == "nl":
            self.branches["BR_SWITCH"] = 2
        elif switching == "bigm":
            self.branches["BR_SWITCH"] = 3
        elif switching == "df":
            # Do nothing, use existing BR_SWITCH values
            pass
        else:
            raise ValueError(
                f"Unknown switching value: {switching}. "
                "Allowed values are 'off' (all lines connected), "
                "'nl' (binary variables, non-linear formulation), "
                "'bigm' (binary variables, Big-M formulation), "
                "'df' (use existing BR_SWITCH), or a numpy array."
            )

    def create_model(self, opf_type="dc", connectivity="off"):
        """Compute the feasible region for the power system.
        Parameters:
        opf_type (str): Type of optimal power flow ('dc', 'acrect', 'acjabr')
        connectivity (str): Connectivity for topology solutions ('off', 'on')
        """
        print(f"=======Computing feasible region ({opf_type}) with connectivity {connectivity}")
        self.ampl = AMPL()
        self.ampl.reset()
        self.ampl.read(Path(__file__).parent / "opf.mod")
        self.ampl.set_data(self.buses, "N")
        self.ampl.set_data(self.generators, "G")
        self.ampl.set_data(self.branches, "L")
        self.ampl.set_data(self.gencost)
        self.ampl.param["CF"] = array2dict(self.cf)
        self.ampl.param["CT"] = array2dict(self.ct)
        self.ampl.param["CG"] = array2dict(self.cg)
        self.ampl.param["OPF_TYPE"] = opf_type
        self.ampl.param["CONNECTIVITY"] = connectivity
        self.ampl.param["BASEMVA"] = self.baseMVA

    def solve_model(self, solver="gurobi", options=""):
        """Solve the model using the specified solver.
        Parameters:
        solver (str): Solver to use ('gurobi', 'cplex', 'cbc')
        options (str): Options for the solver
        Returns:
        dict: Results of the optimization
        """
        print(f"=======Solving model with solver {solver} and options {options}")
        self.ampl.option[solver + "_options"] = options
        self.ampl.solve(solver=solver)

    def solve_opf(self, opf_type="dc", switching="off", connectivity="off", solver="gurobi", options=""):
        """Solve the optimal power flow problem using AMPL.
        Parameters:
        opf_type (str): Type of optimal power flow ('dc', 'acrect', 'acjabr')
        switching (str): Switching strategy ('off', 'nl', 'bigm')
        connectivity (str): Connectivity for topology solutions ('off', 'on')
        solver (str): Solver to use ('gurobi', 'cplex', 'cbc')
        options (str): Options for the solver
        Returns:
        dict: Results of the optimal power flow problem
        """
        self.set_switching(switching)
        self.create_model(opf_type, connectivity)
        self.ampl.eval("minimize objective: sum {g in G} (COST_2[g] * (BASEMVA*Pg[g])^2 + COST_1[g] * (BASEMVA*Pg[g]) + COST_0[g]);")
        self.ampl.eval("option presolve_eps 1e-10;")
        self.solve_model(solver, options)
        return self.get_results_opf(opf_type)

    def get_results_opf(self, opf_type="dc"):
        """Get results from the solved model.
        Returns:
        object: Results of the optimization with attributes like obj, time, generators, buses, branches, etc.
        """
        solver_status = self.ampl.solve_result

        # Create a simple object to hold results
        results = type("Results", (object,), {})()

        # Get the generation results
        Pg = self.ampl.get_variable("Pg").get_values().to_pandas().values.flatten()
        Qg = self.ampl.get_variable("Qg").get_values().to_pandas().values.flatten()

        # Avoid division by zero for Pg_viol
        pmax_pmin_diff = self.generators["PMAX"].values - self.generators["PMIN"].values
        Pg_viol = np.where(
            pmax_pmin_diff == 0,
            0,
            100 * np.maximum(0, Pg - self.generators["PMAX"].values, self.generators["PMIN"].values - Pg) / pmax_pmin_diff,
        )

        # Avoid division by zero for Qg_viol
        qmax_qmin_diff = self.generators["QMAX"].values - self.generators["QMIN"].values
        Qg_viol = np.where(
            qmax_qmin_diff == 0,
            0,
            100 * np.maximum(0, Qg - self.generators["QMAX"].values, self.generators["QMIN"].values - Qg) / qmax_qmin_diff,
        )

        results.generators = pd.DataFrame(
            {"Pg": Pg, "Qg": Qg, "Pg_viol": Pg_viol, "Qg_viol": Qg_viol},
            index=self.ampl.get_variable("Pg").get_values().to_pandas().index,
        )

        # Get the line results
        status = self.ampl.get_variable("status").get_values().to_pandas().values.flatten()
        Pf = self.ampl.get_variable("Pf").get_values().to_pandas().values.flatten()
        Pfa = self.ampl.get_variable("Pfa").get_values().to_pandas().values.flatten()
        Pt = self.ampl.get_variable("Pt").get_values().to_pandas().values.flatten()
        Pta = self.ampl.get_variable("Pta").get_values().to_pandas().values.flatten()
        Qf = self.ampl.get_variable("Qf").get_values().to_pandas().values.flatten()
        Qfa = self.ampl.get_variable("Qfa").get_values().to_pandas().values.flatten()
        Qt = self.ampl.get_variable("Qt").get_values().to_pandas().values.flatten()
        Qta = self.ampl.get_variable("Qta").get_values().to_pandas().values.flatten()
        Sf = Pf + 1j * Qf
        St = Pt + 1j * Qt
        Sf_viol = (
            100
            * np.maximum(0, abs(Sf) - self.branches["RATE_A"].values, -self.branches["RATE_A"].values - abs(Sf))
            / (2 * self.branches["RATE_A"].values)
        )
        St_viol = (
            100
            * np.maximum(0, abs(St) - self.branches["RATE_A"].values, -self.branches["RATE_A"].values - abs(St))
            / (2 * self.branches["RATE_A"].values)
        )
        results.branches = pd.DataFrame(
            {
                "status": status,
                "Pf": Pf,
                "Pt": Pt,
                "Qf": Qf,
                "Qt": Qt,
                "Sf": abs(Sf),
                "St": abs(St),
                "Sf_viol": Sf_viol,
                "St_viol": St_viol,
                "Pfa": Pfa,
                "Pta": Pta,
                "Qfa": Qfa,
                "Qta": Qta,
            },
            index=self.ampl.get_variable("status").get_values().to_pandas().index,
        )

        # Get the voltage results
        if opf_type == "acrect":
            volr = self.ampl.get_variable("Vr").get_values().to_pandas().values.flatten()
            voli = self.ampl.get_variable("Vi").get_values().to_pandas().values.flatten()
            Vm = np.sqrt(volr**2 + voli**2)
            Va = np.arctan2(voli, volr)
        elif opf_type == "acjabr":
            vol2 = self.ampl.get_variable("V2").get_values().to_pandas().values.flatten()
            Vm = np.sqrt(vol2)
            vfvtcosft = self.ampl.get_variable("cosft").get_values().to_pandas().values.flatten()
            vfvt = np.array([Vm[int(self.branches.loc[i, "F_BUS"])] * Vm[int(self.branches.loc[i, "T_BUS"])] for i in range(self.nlin)])
            cosft = np.maximum(-1, np.minimum(1, vfvtcosft / vfvt))
            # Compute angles for all buses
            Va = np.full(self.nbus, np.nan)  # Initialize angles with NaN
            Va[0] = 0  # Reference bus angle is 0
            # Iteratively compute angles
            visited = {0}  # Start with the reference bus
            while len(visited) < self.nbus:
                for line_index in range(self.nlin):
                    f_bus = int(self.branches.loc[line_index, "F_BUS"])
                    t_bus = int(self.branches.loc[line_index, "T_BUS"])
                    if f_bus in visited and np.isnan(Va[t_bus]):
                        Va[t_bus] = Va[f_bus] + np.arccos(cosft[line_index])
                        visited.add(t_bus)
                    elif t_bus in visited and np.isnan(Va[f_bus]):
                        Va[f_bus] = Va[t_bus] - np.arccos(cosft[line_index])
                        visited.add(f_bus)
        else:
            Vm = self.ampl.get_variable("Vm").get_values().to_pandas().values.flatten()
            Va = self.ampl.get_variable("Va").get_values().to_pandas().values.flatten()
        Vm_viol = (
            100
            * np.maximum(0, Vm - self.buses["VMAX"].values, self.buses["VMIN"].values - Vm)
            / (self.buses["VMAX"].values - self.buses["VMIN"].values)
        )
        Va_viol = (
            100
            * np.maximum(0, Va - self.buses["AMAX"].values, self.buses["AMIN"].values - Va)
            / (self.buses["AMAX"].values - self.buses["AMIN"].values)
        )

        # Computation of power injections
        Sd = self.buses["PD"].values + 1j * self.buses["QD"].values
        Sg = Pg + 1j * Qg
        Ssh = self.buses["GS"].values * Vm**2 - 1j * self.buses["BS"].values * Vm**2
        S_viol = Sg @ self.cg - Sd - Ssh - Sf @ self.cf - St @ self.ct
        P_viol = 100 * np.real(S_viol) / sum(self.buses["PD"].values)
        Q_viol = 100 * np.imag(S_viol) / sum(self.buses["QD"].values)

        results.buses = pd.DataFrame(
            {
                "Vm": Vm,
                "Va": Va,
                "Vm_viol": Vm_viol,
                "Va_viol": Va_viol,
                "P_viol": P_viol,
                "Q_viol": Q_viol,
            },
            index=self.ampl.get_variable("Vm").get_values().to_pandas().index,
        )

        # Set other attributes
        results.obj = self.ampl.get_objective("objective").value()
        try:
            results.bestbound = self.ampl.get_value("objective.bestbound")
        except RuntimeError:
            results.bestbound = None
        results.time = self.ampl.get_value("_solve_time")
        results.total_solve_time = self.ampl.get_value("_total_solve_time")
        results.solve_system_time = self.ampl.get_value("_solve_system_time")
        results.solver_status = solver_status
        results.max_viol = float(
            max(
                np.max(np.abs(Pg_viol)),
                np.max(np.abs(Qg_viol)),
                np.max(np.abs(Vm_viol)),
                np.max(np.abs(Va_viol)),
                np.max(np.abs(Sf_viol)),
                np.max(np.abs(St_viol)),
                np.max(np.abs(P_viol)),
                np.max(np.abs(Q_viol)),
            )
        )

        return results

    def is_feasible(self, voltages, angles):
        """
        Check feasibility of a given voltage and angle vector.
        Returns a dictionary with slacks for each constraint.
        If a variable is feasible, slack is 0. Otherwise, it is the absolute value needed to make it feasible.
        """
        slacks = {}

        # 1. Voltage magnitude and angle bounds
        vmin = self.buses["VMIN"].values
        vmax = self.buses["VMAX"].values
        amin = self.buses["AMIN"].values
        amax = self.buses["AMAX"].values

        Vm_slack = np.maximum(0, vmin - voltages) + np.maximum(0, voltages - vmax)
        Va_slack = np.maximum(0, amin - angles) + np.maximum(0, angles - amax)
        slacks["Vm_slack"] = Vm_slack
        slacks["Va_slack"] = Va_slack

        # 2. Compute bus voltages
        v = voltages * np.exp(1j * angles)

        # 3. Compute branch flows
        sf = (self.cf @ v) * np.conj(self.yf @ v)
        st = (self.ct @ v) * np.conj(self.yt @ v)

        # 4. Branch flow bounds (use RATE_A)
        rate_a = self.branches["RATE_A"].values
        abs_sf = np.abs(sf)
        abs_st = np.abs(st)
        Pf_slack = np.maximum(0, abs_sf - rate_a)
        Pt_slack = np.maximum(0, abs_st - rate_a)
        slacks["Sf_slack"] = Pf_slack
        slacks["St_slack"] = Pt_slack

        # 5. Power balance at buses: Pg + jQg = Sbus = v * conj(yb @ v) + Sload
        sd = self.buses["PD"].values + 1j * self.buses["QD"].values
        sb = v * np.conj(self.yb @ v)
        sg = sb + sd  # total generation at each bus

        # 6. Split generator injections among generators at each bus
        pg = np.zeros(self.ngen)
        qg = np.zeros(self.ngen)
        for bus in range(self.nbus):
            gen_indices = self.generators[self.generators["GEN_BUS"] == bus].index
            if len(gen_indices) > 0:
                pmax_total = self.generators.loc[gen_indices, "PMAX"].sum()
                if pmax_total > 0:
                    pg[gen_indices] = np.real(sg[bus]) * self.generators.loc[gen_indices, "PMAX"] / pmax_total
                    qg[gen_indices] = np.imag(sg[bus]) * self.generators.loc[gen_indices, "PMAX"] / pmax_total

        # 7. Generator bounds
        pmin = self.generators["PMIN"].values
        pmax = self.generators["PMAX"].values
        qmin = self.generators["QMIN"].values
        qmax = self.generators["QMAX"].values

        Pg_slack = np.maximum(0, pmin - pg) + np.maximum(0, pg - pmax)
        Qg_slack = np.maximum(0, qmin - qg) + np.maximum(0, qg - qmax)
        slacks["Pg_slack"] = Pg_slack
        slacks["Qg_slack"] = Qg_slack

        return slacks
