import subprocess


def test_main():
    assert subprocess.check_output(["amplpower", "foo", "foobar"], text=True) == "foobar\n"
