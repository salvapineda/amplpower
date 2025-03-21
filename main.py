import logging

from amplpower import PowerSystem

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s", filename="logfile.log", filemode="w")

for switching in ["off", "bigm", "nl"]:
    if switching == "off":
        ps = PowerSystem("./src/amplpower/data/case9.m")
    else:
        ps = PowerSystem("./src/amplpower/data/case9_switching.m")
    for opf_type in ["dc", "acjabr", "acrect", "acpolar"]:
        for solver in ["gurobi", "xpress", "cplex", "copt", "mosek", "highs", "scip", "cbc", "ipopt"]:
            results = ps.solve_opf(
                opf_type=opf_type, switching=switching, connectivity="off", solver=solver, options="outlev=1 timelimit=10"
            )
            logging.info(
                f"Switching: {switching}, OPF Type: {opf_type}, Solver: {solver}, Status: {results['status']}, Objective: {results['obj']}, Time: {results['time']}"
            )
