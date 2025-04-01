=====
Usage
=====

To use the project, follow these steps:

.. code-block:: python

    # Import the PowerSystem class from the amplpower package
    from amplpower import PowerSystem

    # Load the case file
    ps = PowerSystem("case9.m")

    # Solve different types of Optimal Power Flow (OPF):

    # 1. Solve the DC OPF
    dc_results = ps.solve_opf(opf_type='dc')

    # 2. Solve the AC OPF (polar formulation)
    acpolar_results = ps.solve_opf(opf_type='acpolar')

    # 3. Solve the AC OPF (rectangular formulation)
    acrect_results = ps.solve_opf(opf_type='acrect')

    # 4. Solve the Jabr OPF relaxation
    acjabr_results = ps.solve_opf(opf_type='acjabr')
