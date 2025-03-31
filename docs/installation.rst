============
Installation
============

1. **Install amplpower**
   Run the following command to install `amplpower`::

       python -m pip install amplpower

   The `amplpy` package will be installed automatically as a dependency.

2. **Download solver modules**
   After installing `amplpower`, download the solver modules you plan to use. It is important to install at least one solver module, as this also downloads the base module required for AMPL. Without this, you may encounter errors such as "ampl not found". For example, to install solvers like Gurobi, HiGHS, CBC, and IPOPT, run::

       python -m amplpy.modules install gurobi highs coin

   Replace the solver names with the ones relevant to your use case.

3. **Activate an AMPL license**
   By default, `amplpy` includes a demo license with limitations. You can obtain a Community Edition license from [AMPL's portal](https://portal.ampl.com) and activate it using the following command::

       python -m amplpy.modules activate <license-uuid>

   Replace `<license-uuid>` with the UUID of your license.
