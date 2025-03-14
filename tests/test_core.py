from amplpower import PowerSystem
from amplpower import compute


def test_compute():
    assert compute(["a", "bc", "abc"]) == "abc"


def test_powersystem():
    ps = PowerSystem("case9")
    assert ps is not None
