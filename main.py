from amplpower import PowerSystem

ps = PowerSystem("./src/amplpower/data/case9.m")
results = ps.solve_opf(opf_type="dc", switching="off", solver="gurobi", time_limit=60)
