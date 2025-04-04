Changelog
=========

0.0.16 (2025-04-04)
-------------------

* Modification of default Bigm values for AC OTS.

0.0.15 (2025-04-04)
-------------------

* Bus mapping for some pglib networks
* Compute Bigm for AC OTS only if needed

0.0.14 (2025-04-04)
-------------------

* Bug in solver status

0.0.13 (2025-04-04)
-------------------

* Bug options solver

0.0.12 (2025-04-04)
-------------------

* Modify option solver input

0.0.11 (2025-04-01)
-------------------

* Return results for any status

0.0.10 (2025-04-01)
-------------------

* Add formulations of OPF to docs

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
