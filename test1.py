import random  # Import random for random sampling

import numpy as np
import pandas as pd

from amplpower import PowerSystem


def close_lines(ps, edge_index, level=1):
    """Find the indexes of edges close to the given edge up to a specified level.
    Parameters:
    ps (PowerSystem): An instance of the PowerSystem class.
    edge_index (int): Index of the edge for which to find close edges.
    level (int): Depth of adjacency to consider (1 for adjacent, 2 for adjacent of adjacent, etc.).
    Returns:
    list: Indexes of edges connected to the same f_bus or t_bus up to the specified level, excluding the given edge.
    """
    visited_edges = set()
    current_edges = {edge_index}

    for _ in range(level):
        next_edges = set()
        for edge in current_edges:
            branch = ps.branches.iloc[edge]
            f_bus = int(branch["F_BUS"])
            t_bus = int(branch["T_BUS"])
            close_edges = ps.branches[
                (
                    (ps.branches["F_BUS"] == f_bus)
                    | (ps.branches["T_BUS"] == f_bus)
                    | (ps.branches["F_BUS"] == t_bus)
                    | (ps.branches["T_BUS"] == t_bus)
                )
                & (~ps.branches.index.isin(visited_edges))
            ].index.tolist()
            next_edges.update(close_edges)
        visited_edges.update(current_edges)
        current_edges = next_edges

    visited_edges.discard(edge_index)
    return list(visited_edges)


def update_bigm(ps, opf_type="dc", connectivity="off", solver="gurobi", relax="all"):
    """Update the bounds for each line in the power system using a combined approach.
    Parameters:
    ps (PowerSystem): An instance of the PowerSystem class.
    opf_type (str): Type of optimal power flow ('dc', 'acrect', 'acjabr').
    connectivity (str): Connectivity for topology solutions ('off', 'on').
    solver (str): Solver to use ('gurobi', 'cplex', 'cbc').
    relax (str, int, or None): Relaxation strategy for binary variables:
        - "all": Relax all binary variables (default).
        - Any integer: Relax binary variables for lines not in `close_lines` at the specified level.
        - "randX": Randomly select X% of lines, and relax integrality for all other lines (e.g., "rand10" for 10%).
        - "gurobi": Do not relax any binary variables, but use strict mip options and use bestbound.
        - None: Do not update Big-M bounds.
    """
    if relax is None:
        print("No Big-M update performed as relax is None.")
        return 0  # Return 0 as no time is spent updating bounds.
    else:
        total_time = 0
        options = "outlev=0 threads=1"
        for lin_index in range(ps.nlin):
            # Maximize Pfa
            ps.create_model(opf_type, connectivity)
            ps.ampl.eval(f"fix status[{lin_index}]:=0;")
            ps.ampl.eval(f"maximize newbound: Pfa[{lin_index}];")

            options = "outlev=0 threads=1"
            if relax == "all":
                ps.ampl.eval("option relax_integrality 1;")
            elif isinstance(relax, int) and relax > 0:
                relaxed_lines = close_lines(ps, lin_index, level=relax)
                for lin_index in range(ps.nlin):
                    if lin_index not in relaxed_lines:
                        ps.ampl.eval(f"let status[{lin_index}].relax := 1;")
            elif isinstance(relax, str) and relax.startswith("rand"):
                try:
                    percentage = int(relax[4:])
                    if 1 <= percentage <= 99:
                        num_lines_to_keep = max(1, int((percentage / 100) * ps.nlin))
                        selected_lines = set(random.sample(range(ps.nlin), num_lines_to_keep))
                        for line in range(ps.nlin):
                            if line not in selected_lines:
                                ps.ampl.eval(f"let status[{line}].relax := 1;")
                    else:
                        raise ValueError("Percentage out of range (1-99).")
                except ValueError:
                    print(f"Invalid relax value: {relax}. Skipping relaxation.")
                    continue
            elif relax == "gurobi":
                options = "outlev=0 threads=1 timelimit=1 mip:bestbound=1"

            ps.solve_model(solver, options)
            solve_status = ps.ampl.solve_result
            total_time += ps.ampl.get_value("_solve_time")
            if solve_status == "solved":
                if relax == "gurobi":
                    new_bound = ps.ampl.get_value("newbound.bestbound")
                else:
                    new_bound = ps.ampl.get_objective("newbound").value()
                current_bound = ps.branches.loc[lin_index, "PFUPDC"]
                print(f"Line {lin_index} UP: {current_bound} -> {new_bound}")
                ps.branches.loc[lin_index, "PFUPDC"] = np.ceil(100 * new_bound) / 100
            elif solve_status == "infeasible":
                print(f"Infeasible for line {lin_index}")

            # Minimize Pfa
            ps.create_model(opf_type, connectivity)
            ps.ampl.eval(f"fix status[{lin_index}]:=0;")
            ps.ampl.eval(f"minimize newbound: Pfa[{lin_index}];")

            options = "outlev=0 threads=1"
            if relax == "all":
                ps.ampl.eval("option relax_integrality 1;")
            elif isinstance(relax, int) and relax > 0:
                relaxed_lines = close_lines(ps, lin_index, level=relax)
                for lin_index in range(ps.nlin):
                    if lin_index not in relaxed_lines:
                        ps.ampl.eval(f"let status[{lin_index}].relax := 1;")
            elif isinstance(relax, str) and relax.startswith("rand"):
                try:
                    percentage = int(relax[4:])
                    if 1 <= percentage <= 99:
                        num_lines_to_keep = max(1, int((percentage / 100) * ps.nlin))
                        selected_lines = set(random.sample(range(ps.nlin), num_lines_to_keep))
                        for line in range(ps.nlin):
                            if line not in selected_lines:
                                ps.ampl.eval(f"let status[{line}].relax := 1;")
                    else:
                        raise ValueError("Percentage out of range (1-99).")
                except ValueError:
                    print(f"Invalid relax value: {relax}. Skipping relaxation.")
                    continue
            elif relax == "gurobi":
                options = "outlev=0 threads=1 timelimit=1 mip:bestbound=1"

            ps.solve_model(solver, options)
            solve_status = ps.ampl.solve_result
            total_time += ps.ampl.get_value("_solve_time")
            if solve_status == "solved":
                if relax == "gurobi":
                    new_bound = ps.ampl.get_value("newbound.bestbound")
                else:
                    new_bound = ps.ampl.get_objective("newbound").value()
                current_bound = ps.branches.loc[lin_index, "PFLODC"]
                print(f"Line {lin_index} LO: {current_bound} -> {new_bound}")
                ps.branches.loc[lin_index, "PFLODC"] = np.floor(100 * new_bound) / 100
            elif solve_status == "infeasible":
                print(f"Infeasible for line {lin_index}")

        print(f"Total time for updating bounds: {total_time:.2f} seconds")
        return total_time


# Iterate over relax values
results_df = pd.DataFrame(columns=["relax", "time_update_bigm", "objective", "ots_time"])
# for relax in [ "all", 1, 2, "rand10", "rand20"]:
for relax in ["gurobi", 1, 2, "rand10", "rand20"]:
    print(f"Running update_bigm with relax={relax}")

    # Create a new PowerSystem instance for each iteration
    ps = PowerSystem("./src/amplpower/data/case118Blumsack.m")
    ps.set_switching("bigm")
    ps.ubcost = 1600

    # Update bounds and measure time
    time_update_bigm = update_bigm(ps, relax=relax)

    # Solve OTS and compute results
    ps.compute_ub_cost()
    results_ots = ps.solve_opf(opf_type="dc", switching="bigm", solver="gurobi", options="outlev=1 threads=1 timelimit=100")

    # Append results to the DataFrame
    results_df = pd.concat(
        [
            results_df,
            pd.DataFrame(
                {"relax": [relax], "time_update_bigm": [time_update_bigm], "objective": [results_ots.obj], "ots_time": [results_ots.time]}
            ),
        ],
        ignore_index=True,
    )

# Display results on screen
print(results_df)
