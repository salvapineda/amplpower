import sys  # Import the sys module
from pathlib import Path

import numpy as np

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


def new_bound(ps, opf_type, connectivity, solver, direction, lin_index, fix_status, relaxed_lines):
    """
    Solves for the maximum or minimum value of Pfa[lin_index] with specified relaxed lines.

    Parameters:
    ps (PowerSystem): Power system instance.
    opf_type (str): OPF type.
    connectivity (str): Connectivity option.
    solver (str): Solver name.
    direction (str): "max" for maximization, "min" for minimization.
    lin_index (int): Index of the line to optimize.
    fix_status (int): 1 to force the line to be connected, 0 to force it to be disconnected.
    relaxed_lines (list): List of line indices whose binary variables are relaxed.

    Returns:
    tuple: (new_bound_value, solve_status, solve_time)
    """
    ps.create_model(opf_type, connectivity)
    ps.ampl.eval(f"fix status[{lin_index}]:={fix_status};")
    if direction == "max":
        ps.ampl.eval(f"maximize newbound: Pfa[{lin_index}];")
    else:
        ps.ampl.eval(f"minimize newbound: Pfa[{lin_index}];")
    ps.ampl.eval(
        f"subject to upper_bound_cost: sum {{g in G}} (COST_2[g] * (BASEMVA*Pg[g])^2 + COST_1[g] * (BASEMVA*Pg[g]) + COST_0[g]) <= {ps.ubcost};"
    )
    # Relax integrality for specified lines
    for line in range(ps.nlin):
        if line in relaxed_lines:
            ps.ampl.eval(f"let status[{line}].relax := 1;")
    options = "outlev=0 threads=1"
    ps.solve_model(solver, options)
    solve_status = ps.ampl.solve_result
    solve_time = ps.ampl.get_value("_solve_time")
    if solve_status == "solved":
        new_bound_value = ps.ampl.get_objective("newbound").value()
    else:
        new_bound_value = None
    return new_bound_value, solve_status, solve_time


def update_bigm(ps, opf_type="dc", connectivity="off", solver="gurobi", sorted="default", relax="all"):
    """
    Update the bounds for each line in the power system using a combined approach.

    Parameters:
    ps (PowerSystem): An instance of the PowerSystem class.
    opf_type (str): Type of optimal power flow ('dc', 'acrect', 'acjabr').
    connectivity (str): Connectivity for topology solutions ('off', 'on').
    solver (str): Solver to use ('gurobi', 'cplex', 'cbc').
    sorted (str): Determines line order:
        - "default": use range(ps.nlin)
        - "ranking": use ranking_lines(ps)
        - "random": use rand_line(ps) ps.nlin times (with replacement)
    relax (str): Strategy for relaxing binary variables:
        - "all": Relax all binary variables except the current line.
        - "close1": Relax binary variables for lines not in close_lines at level 1.
        - "close2": Relax binary variables for lines not in close_lines at level 2.
    Returns:
    tuple: (total_time, avg_improvement)
    """
    total_time = 0
    improvements = []
    updated_bounds = []  # Store (line, type, f_bus, t_bus, old, new, improvement)
    # Determine line order
    if sorted == "ranking":
        line_order = ranking_lines(ps)
    elif sorted == "random":
        line_order = [rand_line(ps) for _ in range(ps.nlin)]
    else:  # "default"
        line_order = range(ps.nlin)
    # First, update PFMAX and PFMIN with fix_status=1 for all lines
    for lin_index in line_order:
        print(f"Processing line {lin_index} (sorted={sorted}, relax={relax})")
        # Determine relaxed lines based on relax parameter
        if relax == "all":
            relaxed_lines = [i for i in range(ps.nlin) if i != lin_index]
        elif relax.startswith("close"):
            try:
                level = int(relax.replace("close", ""))
            except Exception:
                print(f"Invalid relax: {relax}. Skipping line {lin_index}.")
                continue
            close = set(close_lines(ps, lin_index, level=level))
            relaxed_lines = [i for i in range(ps.nlin) if i not in close]
        else:
            print(f"Unknown relax: {relax}. Skipping line {lin_index}.")
            continue

        f_bus = ps.branches.loc[lin_index, "F_BUS"]
        t_bus = ps.branches.loc[lin_index, "T_BUS"]

        # Maximize Pfa with fix_status=1 (connected, update PFMAX)
        max_bound_1, max_status_1, max_time_1 = new_bound(ps, opf_type, connectivity, solver, "max", lin_index, 1, relaxed_lines)
        total_time += max_time_1
        current_bound_up_1 = ps.branches.loc[lin_index, "PFMAX"]
        if max_status_1 == "solved" and max_bound_1 is not None and current_bound_up_1 != 0:
            rounded_max_bound_1 = np.ceil(1e10 * max_bound_1) / 1e10
            ps.branches.loc[lin_index, "PFMAX"] = rounded_max_bound_1
            improvement_up_1 = abs((rounded_max_bound_1 - current_bound_up_1) / current_bound_up_1) * 100
            print(f"Line {lin_index} UP (fix_status=1): {current_bound_up_1} -> {rounded_max_bound_1} ({improvement_up_1:.2f}%)")
            if improvement_up_1 > 0.01:
                updated_bounds.append(("UP-1", lin_index, f_bus, t_bus, current_bound_up_1, rounded_max_bound_1, improvement_up_1))
        else:
            improvement_up_1 = 0.0
            if max_status_1 == "infeasible":
                print(f"Infeasible for line {lin_index} (max, fix_status=1)")
        improvements.append(improvement_up_1)

        # Minimize Pfa with fix_status=1 (connected, update PFMIN)
        min_bound_1, min_status_1, min_time_1 = new_bound(ps, opf_type, connectivity, solver, "min", lin_index, 1, relaxed_lines)
        total_time += min_time_1
        current_bound_lo_1 = ps.branches.loc[lin_index, "PFMIN"]
        if min_status_1 == "solved" and min_bound_1 is not None and current_bound_lo_1 != 0:
            rounded_min_bound_1 = np.floor(1e10 * min_bound_1) / 1e10
            ps.branches.loc[lin_index, "PFMIN"] = rounded_min_bound_1
            improvement_lo_1 = abs((rounded_min_bound_1 - current_bound_lo_1) / current_bound_lo_1) * 100
            print(f"Line {lin_index} LO (fix_status=1): {current_bound_lo_1} -> {rounded_min_bound_1} ({improvement_lo_1:.2f}%)")
            if improvement_lo_1 > 0.01:
                updated_bounds.append(("LO-1", lin_index, f_bus, t_bus, current_bound_lo_1, rounded_min_bound_1, improvement_lo_1))
        else:
            improvement_lo_1 = 0.0
            if min_status_1 == "infeasible":
                print(f"Infeasible for line {lin_index} (min, fix_status=1)")
        improvements.append(improvement_lo_1)

    # Then, update PFUPDC and PFLODC with fix_status=0 for all lines
    for lin_index in line_order:
        print(f"Processing line {lin_index} (sorted={sorted}, relax={relax})")
        if relax == "all":
            relaxed_lines = [i for i in range(ps.nlin) if i != lin_index]
        elif relax.startswith("close"):
            try:
                level = int(relax.replace("close", ""))
            except Exception:
                print(f"Invalid relax: {relax}. Skipping line {lin_index}.")
                continue
            close = set(close_lines(ps, lin_index, level=level))
            relaxed_lines = [i for i in range(ps.nlin) if i not in close]
        else:
            print(f"Unknown relax: {relax}. Skipping line {lin_index}.")
            continue

        f_bus = ps.branches.loc[lin_index, "F_BUS"]
        t_bus = ps.branches.loc[lin_index, "T_BUS"]

        # Maximize Pfa with fix_status=0 (disconnected, update PFUPDC)
        max_bound_0, max_status_0, max_time_0 = new_bound(ps, opf_type, connectivity, solver, "max", lin_index, 0, relaxed_lines)
        total_time += max_time_0
        current_bound_up = ps.branches.loc[lin_index, "PFUPDC"]
        if max_status_0 == "solved" and max_bound_0 is not None and current_bound_up != 0:
            rounded_max_bound_0 = np.ceil(1e10 * max_bound_0) / 1e10
            ps.branches.loc[lin_index, "PFUPDC"] = rounded_max_bound_0
            improvement_up_0 = abs((rounded_max_bound_0 - current_bound_up) / current_bound_up) * 100
            print(f"Line {lin_index} UP (fix_status=0): {current_bound_up} -> {rounded_max_bound_0} ({improvement_up_0:.2f}%)")
            if improvement_up_0 > 0.01:
                updated_bounds.append(("UP-0", lin_index, f_bus, t_bus, current_bound_up, rounded_max_bound_0, improvement_up_0))
        else:
            improvement_up_0 = 0.0
            if max_status_0 == "infeasible":
                print(f"Infeasible for line {lin_index} (max, fix_status=0)")
        improvements.append(improvement_up_0)

        # Minimize Pfa with fix_status=0 (disconnected, update PFLODC)
        min_bound_0, min_status_0, min_time_0 = new_bound(ps, opf_type, connectivity, solver, "min", lin_index, 0, relaxed_lines)
        total_time += min_time_0
        current_bound_lo = ps.branches.loc[lin_index, "PFLODC"]
        if min_status_0 == "solved" and min_bound_0 is not None and current_bound_lo != 0:
            rounded_min_bound_0 = np.floor(1e10 * min_bound_0) / 1e10
            ps.branches.loc[lin_index, "PFLODC"] = rounded_min_bound_0
            improvement_lo_0 = abs((rounded_min_bound_0 - current_bound_lo) / current_bound_lo) * 100
            print(f"Line {lin_index} LO (fix_status=0): {current_bound_lo} -> {rounded_min_bound_0} ({improvement_lo_0:.2f}%)")
            if improvement_lo_0 > 0.01:
                updated_bounds.append(("LO-0", lin_index, f_bus, t_bus, current_bound_lo, rounded_min_bound_0, improvement_lo_0))
        else:
            improvement_lo_0 = 0.0
            if min_status_0 == "infeasible":
                print(f"Infeasible for line {lin_index} (min, fix_status=0)")
        improvements.append(improvement_lo_0)

    avg_improvement = np.mean(improvements) if improvements else 0.0
    print(f"Total time for updating bounds: {total_time:.2f} seconds")
    print(f"Average percentage improvement: {avg_improvement:.2f}%")
    # Print all updated bounds with improvement > 0%
    print("Updated bounds with improvement > 0%:")
    for bound_type, lin_index, f_bus, t_bus, old, new, imp in updated_bounds:
        print(f"Line {lin_index} ({f_bus}-{t_bus}) {bound_type}: {old} -> {new} ({imp:.2f}%)")
    return total_time, avg_improvement


def ranking_nodes(ps):
    """
    Returns a list of nodes (buses) sorted by the number of lines connected to each node, in descending order.
    Nodes with the most connections appear first, while leaf nodes (with the fewest connections) appear last.

    Parameters:
    ps (PowerSystem): An instance of the PowerSystem class, which must have a 'branches' DataFrame
                      with columns 'F_BUS' and 'T_BUS' representing the from and to buses of each line.

    Returns:
    list: List of node identifiers (bus numbers), sorted from most to least connected.
    """
    # Count connections per node
    f_counts = ps.branches["F_BUS"].value_counts()
    t_counts = ps.branches["T_BUS"].value_counts()
    total_counts = f_counts.add(t_counts, fill_value=0)
    # Sort nodes by number of connections (descending)
    sorted_nodes = total_counts.sort_values(ascending=False).index.tolist()
    return sorted_nodes


def lines_in_node(ps, node):
    """
    Returns a list of indices of lines connected to the specified node.

    Parameters:
    ps (PowerSystem): An instance of the PowerSystem class, which must have a 'branches' DataFrame
                      with columns 'F_BUS' and 'T_BUS' representing the from and to buses of each line.
    node (int): The node (bus) number for which to find all connected lines.

    Returns:
    list: List of indices (integers) of lines connected to the given node.
    """
    mask = (ps.branches["F_BUS"] == node) | (ps.branches["T_BUS"] == node)
    return ps.branches[mask].index.tolist()


def ranking_lines(ps):
    """
    Returns a list of line indices sorted by the importance of the nodes they are connected to.
    Lines connected to higher-ranked nodes (from ranking_nodes) appear first.

    Parameters:
    ps (PowerSystem): An instance of the PowerSystem class.

    Returns:
    list: List of line indices, sorted by node importance.
    """
    node_ranking = ranking_nodes(ps)
    node_rank_map = {node: rank for rank, node in enumerate(node_ranking)}
    line_scores = []
    for idx, branch in ps.branches.iterrows():
        f_rank = node_rank_map.get(branch["F_BUS"], len(node_ranking))
        t_rank = node_rank_map.get(branch["T_BUS"], len(node_ranking))
        # Lower rank value means higher importance, so use min
        line_score = min(f_rank, t_rank)
        line_scores.append((idx, line_score))
    # Sort lines by score (ascending: most important nodes first)
    sorted_lines = [idx for idx, _ in sorted(line_scores, key=lambda x: x[1])]
    return sorted_lines


def rand_line(ps):
    """
    Assigns each line a score equal to the sum of lines connected to its two buses,
    and returns a random line index with probability proportional to its score.

    Parameters:
    ps (PowerSystem): An instance of the PowerSystem class.

    Returns:
    int: Index of the randomly selected line.
    """
    # Count number of lines connected to each bus
    f_counts = ps.branches["F_BUS"].value_counts()
    t_counts = ps.branches["T_BUS"].value_counts()
    total_counts = f_counts.add(t_counts, fill_value=0)

    # Compute points for each line
    points = []
    for _idx, branch in ps.branches.iterrows():
        f_bus = branch["F_BUS"]
        t_bus = branch["T_BUS"]
        score = total_counts.get(f_bus, 0) + total_counts.get(t_bus, 0)
        points.append(score)
    points = np.array(points)

    # Normalize to probabilities
    prob = points / points.sum()
    # Randomly select a line index
    selected_idx = np.random.choice(ps.branches.index, p=prob)
    return selected_idx


def ots_heuristic(ps, opf_type="dc", connectivity="off", solver="gurobi"):
    """
    Runs a simple OTS heuristic: creates the model, sets the objective to minimize total generation cost,
    solves the model, and returns the objective value.

    Parameters:
    ps (PowerSystem): Power system instance.
    opf_type (str): OPF type.
    connectivity (str): Connectivity option.
    solver (str): Solver to use.

    Returns:
    float: Objective value (minimum cost) or None if not solved.
    """
    ps.create_model(opf_type, connectivity)
    ps.ampl.eval("minimize cost: sum {g in G} (COST_2[g] * (BASEMVA*Pg[g])^2 + COST_1[g] * (BASEMVA*Pg[g]) + COST_0[g]);")
    options = "outlev=1 threads=1 timelimit=10 mip:heurfrac=1"
    ps.solve_model(solver, options)
    solve_status = ps.ampl.solve_result
    if solve_status == "solved":
        obj_value = ps.ampl.get_objective("cost").value()
    else:
        obj_value = None
    return obj_value


def random_demand(ps, seed=None):
    """
    Modifies the 'PD' column of ps.buses by multiplying each value by a random number
    drawn from a uniform distribution between 0.9 and 1.1.

    Parameters:
    ps (PowerSystem): Power system instance.
    seed (int): Seed for the random number generator, for reproducibility.
    """
    if seed is not None:
        np.random.seed(seed)
    for b in ps.buses.index:
        ps.buses.loc[b, "PD"] *= np.random.uniform(0.9, 1.1)


def run_ots(m_file, demand_seed, sorted_val, relax_val):
    """
    Runs the OTS solving process with specified parameters.

    Parameters:
    m_file (str): Path to the .m file.
    demand_seed (int): Seed for random demand.
    sorted_val (str): Value for the 'sorted' parameter in update_bigm.
    relax_val (str): Value for the 'relax' parameter in update_bigm.

    Returns:
    tuple: (time_update_bigm, improvement_bigm, ots_time, objective)
    """
    # Create a new PowerSystem instance
    ps = PowerSystem(m_file)
    random_demand(ps, seed=demand_seed)
    ps.set_switching("bigm")
    ps.ubcost = ots_heuristic(ps, opf_type="dc", connectivity="off", solver="gurobi")

    # Update bounds and measure time
    time_update_bigm, improvement_bigm = update_bigm(ps, sorted=sorted_val, relax=relax_val)

    # Solve OTS and compute results
    results_ots = ps.solve_opf(opf_type="dc", switching="bigm", solver="gurobi", options="outlev=1 threads=1 timelimit=1000")

    return time_update_bigm, improvement_bigm, results_ots.time, results_ots.obj


# Define cases
cases = []
for seed in range(100):
    for sorted_val in ["default", "ranking"]:
        for relax_val in ["all", "close1", "close2"]:
            cases.append({"seed": seed, "sorted": sorted_val, "relax": relax_val})

# Define SLURM index
slurm_index = 0  # Example value, you can change this

# Get the case corresponding to the SLURM index
if 0 <= slurm_index < len(cases):
    case = cases[slurm_index]
    seed, sorted_val, relax_val = case["seed"], case["sorted"], case["relax"]
    print(f"Running update_bigm with seed={seed}, sorted={sorted_val}, relax={relax_val}")

    # Run OTS and store results
    time_update_bigm, improvement_bigm, ots_time, objective = run_ots("./src/amplpower/data/case118Blumsack.m", seed, sorted_val, relax_val)

    # Write results to CSV file
    filename = f"{slurm_index}.csv"
    filepath = Path(filename)
    with filepath.open("w") as f:
        f.write("seed,sorted,relax,time_update_bigm,improvement_bigm,ots_time,objective\n")
        f.write(f"{seed},{sorted_val},{relax_val},{time_update_bigm},{improvement_bigm},{ots_time},{objective}\n")

    print(f"Results written to {filename}")
else:
    print(f"SLURM index {slurm_index} is out of range (0-{len(cases) - 1}).")
    sys.exit(1)
