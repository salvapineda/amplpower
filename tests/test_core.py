from amplpower import PowerSystem
from amplpower import compute


def test_compute():
    assert compute(["a", "bc", "abc"]) == "abc"


def test_powersystem():
    ps = PowerSystem("./src/amplpower/data/case9.m")
    ps.solve_opf()
    assert ps is not None
