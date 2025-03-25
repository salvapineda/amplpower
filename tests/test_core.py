import pytest

from amplpower import PowerSystem
from amplpower import compute


def test_compute():
    assert compute(["a", "bc", "abc"]) == "abc"


@pytest.mark.parametrize(
    "opf_type",
    [
        "dc",
        "acjabr",
        "acrect",
        "acpolar",
    ],
)
def test_opf(opf_type):
    ps = PowerSystem("./src/amplpower/data/case9.m")
    results = ps.solve_opf(opf_type=opf_type, switching="off", connectivity="off", solver="gurobi", options="outlev=1 timelimit=5")
    print(results)
    assert results["status"] == "solved"


@pytest.mark.parametrize("opf_type", ["dc", "acjabr", "acrect", "acpolar"])
@pytest.mark.parametrize("switching", ["bigm", "nl"])
@pytest.mark.parametrize("connectivity", ["on", "off"])
def test_ots(opf_type, switching, connectivity):
    ps = PowerSystem("./src/amplpower/data/case9_switching.m")
    results = ps.solve_opf(
        opf_type=opf_type, switching=switching, connectivity=connectivity, solver="gurobi", options="outlev=1 timelimit=5"
    )
    print(results)
    assert results["status"] == "solved"
