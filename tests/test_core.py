import pytest

from amplpower import PowerSystem
from amplpower import compute


def test_compute():
    assert compute(["a", "bc", "abc"]) == "abc"


@pytest.mark.parametrize("solver", ["gurobi", "xpress", "cplex", "copt", "mosek", "highs", "scip", "cbc", "ipopt"])
def test_opf_dc(solver):
    ps = PowerSystem("./src/amplpower/data/case9.m")
    results = ps.solve_opf(opf_type="dc", switching="off", connectivity="off", solver=solver, options="outlev=1 timelimit=5")
    assert results["status"] == "solved"


@pytest.mark.parametrize("solver", ["gurobi", "xpress", "cplex", "mosek", "scip", "ipopt"])
def test_opf_jabr(solver):
    ps = PowerSystem("./src/amplpower/data/case9.m")
    results = ps.solve_opf(opf_type="acjabr", switching="off", connectivity="off", solver=solver, options="outlev=1 timelimit=5")
    assert results["status"] == "solved"


@pytest.mark.parametrize("solver", ["gurobi", "scip", "ipopt"])
def test_opf_rect(solver):
    ps = PowerSystem("./src/amplpower/data/case9.m")
    results = ps.solve_opf(opf_type="acrect", switching="off", connectivity="off", solver=solver, options="outlev=1 timelimit=5")
    assert results["status"] == "solved"


@pytest.mark.parametrize("solver", ["gurobi", "ipopt"])
def test_opf_polar(solver):
    ps = PowerSystem("./src/amplpower/data/case9.m")
    results = ps.solve_opf(opf_type="acpolar", switching="off", connectivity="off", solver=solver, options="outlev=1 timelimit=5")
    assert results["status"] == "solved"


@pytest.mark.parametrize("solver", ["gurobi", "xpress", "cplex", "copt", "mosek", "scip"])
@pytest.mark.parametrize("connectivity", ["on", "off"])
def test_ots_dc(solver, connectivity):
    ps = PowerSystem("./src/amplpower/data/case9_switching.m")
    results = ps.solve_opf(opf_type="dc", switching="bigm", connectivity=connectivity, solver=solver, options="outlev=1 timelimit=5")
    assert results["status"] == "solved"


@pytest.mark.parametrize("solver", ["gurobi", "xpress", "cplex", "copt", "mosek", "scip"])
@pytest.mark.parametrize("connectivity", ["on", "off"])
def test_ots_jabr(solver, connectivity):
    ps = PowerSystem("./src/amplpower/data/case9_switching.m")
    results = ps.solve_opf(opf_type="acjabr", switching="bigm", connectivity=connectivity, solver=solver, options="outlev=1 timelimit=5")
    assert results["status"] == "solved"


@pytest.mark.parametrize("solver", ["gurobi", "scip"])
@pytest.mark.parametrize("connectivity", ["on", "off"])
def test_ots_rect(solver, connectivity):
    ps = PowerSystem("./src/amplpower/data/case9_switching.m")
    results = ps.solve_opf(opf_type="acrect", switching="bigm", connectivity=connectivity, solver=solver, options="outlev=1 timelimit=5")
    assert results["status"] == "solved"


@pytest.mark.parametrize("solver", ["gurobi", "xpress", "cplex", "copt", "mosek", "scip", "ipopt"])
@pytest.mark.parametrize("connectivity", ["on", "off"])
def test_ots_polar(solver, connectivity):
    ps = PowerSystem("./src/amplpower/data/case9_switching.m")
    results = ps.solve_opf(opf_type="acpolar", switching="bigm", connectivity=connectivity, solver=solver, options="outlev=1 timelimit=5")
    assert results["status"] == "solved"


@pytest.mark.parametrize("solver", ["gurobi", "mosek", "scip"])
@pytest.mark.parametrize("connectivity", ["on", "off"])
def test_otsn_dc(solver, connectivity):
    ps = PowerSystem("./src/amplpower/data/case9_switching.m")
    results = ps.solve_opf(opf_type="dc", switching="nl", connectivity=connectivity, solver=solver, options="outlev=1 timelimit=5")
    assert results["status"] == "solved"


@pytest.mark.parametrize("solver", ["gurobi", "scip"])
@pytest.mark.parametrize("connectivity", ["on", "off"])
def test_otsn_jabr(solver, connectivity):
    ps = PowerSystem("./src/amplpower/data/case9_switching.m")
    results = ps.solve_opf(opf_type="acjabr", switching="nl", connectivity=connectivity, solver=solver, options="outlev=1 timelimit=5")
    assert results["status"] == "solved"


@pytest.mark.parametrize("solver", ["gurobi", "scip"])
@pytest.mark.parametrize("connectivity", ["on", "off"])
def test_otsn_rect(solver, connectivity):
    ps = PowerSystem("./src/amplpower/data/case9_switching.m")
    results = ps.solve_opf(opf_type="acrect", switching="nl", connectivity=connectivity, solver=solver, options="outlev=1 timelimit=5")
    assert results["status"] == "solved"


@pytest.mark.parametrize("solver", ["gurobi"])
@pytest.mark.parametrize("connectivity", ["on", "off"])
def test_otsn_polar(solver, connectivity):
    ps = PowerSystem("./src/amplpower/data/case9_switching.m")
    results = ps.solve_opf(opf_type="acpolar", switching="nl", connectivity=connectivity, solver=solver, options="outlev=1 timelimit=5")
    assert results["status"] == "solved"
