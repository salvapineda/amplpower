import pytest

from amplpower import PowerSystem
from amplpower import compute


def test_compute():
    assert compute(["a", "bc", "abc"]) == "abc"


@pytest.mark.parametrize(
    ("opf_type", "solver"),
    [
        ("dc", "highs"),
        ("acjabr", "ipopt"),
        ("acrect", "ipopt"),
        ("acpolar", "ipopt"),
    ],
)
def test_opf(opf_type, solver):
    ps = PowerSystem("./src/amplpower/data/case9.m")
    results = ps.solve_opf(opf_type=opf_type, switching="off", connectivity="off", solver=solver)
    print(f"Objective function value: {results.obj}")
    print(f"Time: {results.time}")
    print(results.generators)
    print(results.buses)
    print(results.branches)
    assert results.solver_status is not None


@pytest.mark.parametrize("opf_type", ["dc", "acjabr", "acrect", "acpolar"])
@pytest.mark.parametrize("switching", ["bigm", "nl"])
def test_ots(opf_type, switching):
    ps = PowerSystem("./src/amplpower/data/case9_switching.m")
    results = ps.solve_opf(opf_type=opf_type, switching=switching, connectivity="off", solver="scip", options="outlev=1 timelimit=5")
    print(f"Objective function value: {results.obj}")
    print(f"Time: {results.time}")
    print(results.generators)
    print(results.buses)
    print(results.branches)
    assert results.solver_status is not None
