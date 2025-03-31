Changelog
=========

0.0.9 (2025-03-26)
------------------

* Tightest computation of bigMs for AC OTS.
* Tightest bounds for variables cosft and sinft.
* Add violations of AC constraints to results.
* Test use open-source solvers.

0.0.8 (2025-03-25)
------------------

* Corrected voltage results for AC jabr relaxation (acjabr).

0.0.7 (2025-03-25)
------------------

* Included default COST2 for generators (gencos) if not provided.
* Corrected voltage results for AC rectangular (acrect).

0.0.6 (2025-03-21)
------------------

* Added support for solving optimal power flow (OPF) problems: DC OPF, AC OPF (both rectangular and polar coordinates) and AC relaxation proposed by Jabr.
* Added functionality for solving the optimal transmission switching (OTS) problem: Big-M and non-linear formulations. Option to include or exclude connectivity constraints for the OTS problem

0.0.0 (2025-03-14)
------------------

* First release on PyPI.
