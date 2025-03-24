from amplpower import PowerSystem
from amplpower import compute


def test_compute():
    assert compute(["a", "bc", "abc"]) == "abc"


def test_opf_dc():
    ps = PowerSystem("./src/amplpower/data/case9.m")
    results = ps.solve_opf(opf_type="dc", switching="off", connectivity="off", solver="gurobi", options="outlev=1 timelimit=5")
    assert results["status"] == "solved"


def test_opf_jabr():
    ps = PowerSystem("./src/amplpower/data/case9.m")
    results = ps.solve_opf(opf_type="acjabr", switching="off", connectivity="off", solver="gurobi", options="outlev=1 timelimit=5")
    assert results["status"] == "solved"


def test_opf_rect():
    ps = PowerSystem("./src/amplpower/data/case9.m")
    results = ps.solve_opf(opf_type="acrect", switching="off", connectivity="off", solver="gurobi", options="outlev=1 timelimit=5")
    assert results["status"] == "solved"


def test_opf_polar():
    ps = PowerSystem("./src/amplpower/data/case9.m")
    results = ps.solve_opf(opf_type="acpolar", switching="off", connectivity="off", solver="gurobi", options="outlev=1 timelimit=5")
    assert results["status"] == "solved"


def test_ots_dc():
    ps = PowerSystem("./src/amplpower/data/case9_switching.m")
    results = ps.solve_opf(opf_type="dc", switching="bigm", connectivity="off", solver="gurobi", options="outlev=1 timelimit=5")
    assert results["status"] == "solved"


def test_ots_jabr():
    ps = PowerSystem("./src/amplpower/data/case9_switching.m")
    results = ps.solve_opf(opf_type="acjabr", switching="bigm", connectivity="off", solver="gurobi", options="outlev=1 timelimit=5")
    assert results["status"] == "solved"


def test_ots_rect():
    ps = PowerSystem("./src/amplpower/data/case9_switching.m")
    results = ps.solve_opf(opf_type="acrect", switching="bigm", connectivity="off", solver="gurobi", options="outlev=1 timelimit=5")
    assert results["status"] == "solved"


def test_ots_polar():
    ps = PowerSystem("./src/amplpower/data/case9_switching.m")
    results = ps.solve_opf(opf_type="acpolar", switching="bigm", connectivity="off", solver="gurobi", options="outlev=1 timelimit=5")
    assert results["status"] == "solved"


def test_otsn_dc():
    ps = PowerSystem("./src/amplpower/data/case9_switching.m")
    results = ps.solve_opf(opf_type="dc", switching="nl", connectivity="off", solver="gurobi", options="outlev=1 timelimit=5")
    assert results["status"] == "solved"


def test_otsn_jabr():
    ps = PowerSystem("./src/amplpower/data/case9_switching.m")
    results = ps.solve_opf(opf_type="acjabr", switching="nl", connectivity="off", solver="gurobi", options="outlev=1 timelimit=5")
    assert results["status"] == "solved"


def test_otsn_rect():
    ps = PowerSystem("./src/amplpower/data/case9_switching.m")
    results = ps.solve_opf(opf_type="acrect", switching="nl", connectivity="off", solver="gurobi", options="outlev=1 timelimit=5")
    assert results["status"] == "solved"


def test_otsn_polar():
    ps = PowerSystem("./src/amplpower/data/case9_switching.m")
    results = ps.solve_opf(opf_type="acpolar", switching="nl", connectivity="off", solver="gurobi", options="outlev=1 timelimit=5")
    assert results["status"] == "solved"
