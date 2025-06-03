Changelog
=========

0.0.31 (2025-06-03)
-------------------

* Remove minus infinity from variables in AMPL model

0.0.30 (2025-06-03)
-------------------

* Floor and ceil function for initial BigM values in OTS

0.0.29 (2025-06-02)
-------------------

* Bug with bound for power flows in OTS problem

0.0.28 (2025-05-22)
-------------------

* Upper and lower bound for power flows

0.0.27 (2025-05-20)
-------------------

* Remove cost definition from the model due to bad scaling in constraints

0.0.26 (2025-05-19)
-------------------

* Bug with try except in results

0.0.25 (2025-05-07)
-------------------

* Split solve_opf function
* Add upper bound to total generation cost
* Split creat model and solve model for clarity
* Change result output
* The AMPL model is now stored in self.ampl and can be modified by user
* Add best bound to results if available

0.0.24 (2025-05-05)
-------------------

* No default options for solver

0.0.23 (2025-04-24)
-------------------

* Fix initialize generation leven when multiple units at the same bus

0.0.22 (2025-04-24)
-------------------

* Fix division by zero in generator violation calculation
* Chage options input to the solver

0.0.20 (2025-04-11)
-------------------

* Add maximum violation of AC constraints to results.

0.0.19 (2025-04-09)
-------------------

* Bug bus mapping

0.0.18 (2025-04-08)
-------------------

* BigM for AC OTS are computed exploring all critical points
* Compute bounds for real and imaginary parts of voltage

0.0.17 (2025-04-04)
-------------------

* Bug relatex to BR_X negative in some networks.

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
