=====
Usage
=====

To use the project:

.. code-block:: python

    # Import the PowerSystem class from the amplpower package
    from amplpower import PowerSystem

.. code-block:: python

    # Load the case file
    ps = PowerSystem("case9.m")

.. code-block:: python

    # Solve the DC OPF with specified parameters
    results = ps.solve_opf(opf_type='dc', switching='off', connectivity='off', solver='gurobi', options='outlev=1 timelimit=60')
