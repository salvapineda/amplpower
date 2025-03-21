=====
Usage
=====

To use the project:

.. code-block:: python

    # Import the PowerSystem class from the amplpower package
    from amplpower import PowerSystem

    # Load the case file
    ps = PowerSystem("case9.m")

    # Solve the DC OPF
    results = ps.solve_opf(opf_type='dc', switching='off', connectivity='off', solver='gurobi'')

    # Solve the AC OPF (rectangular formulation)
    results = ps.solve_opf(opf_type='acrect', switching='off', connectivity='off', solver='gurobi')

    # Solve the AC OPF (polar formulation)
    results = ps.solve_opf(opf_type='acpolar', switching='off', connectivity='off', solver='gurobi')

    # Solve the Jabr OPF relaxation
    results = ps.solve_opf(opf_type='acjabr', switching='off', connectivity='off', solver='gurobi')

    # Solve the DC OTS using the bigM formulation without imposing connectivity constraints
    results = ps.solve_ots(ots_type='dc', switching='bigm', connectivity='off', solver='gurobi')

    # Solve the DC OTS using the non-linear formulation without imposing connectivity constraints
    results = ps.solve_ots(ots_type='dc', switching='nl', connectivity='off', solver='gurobi')

    # Solve the DC OTS using the bigM formulation and imposing connectivity constraints
    results = ps.solve_ots(ots_type='dc', switching='bigm', connectivity='on', solver='gurobi')

    # Solve the DC OTS using the non-linear formulation and imposing connectivity constraints
    results = ps.solve_ots(ots_type='dc', switching='nl', connectivity='on', solver='gurobi')
