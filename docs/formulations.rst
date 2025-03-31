===========================
Formulations and Models
===========================

This section provides an overview of the formulations and optimization models used in the package.

Nomenclature
------------

Indexes
~~~~~~~
.. list-table::
   :widths: 20 80
   :header-rows: 0

   * - :math:`n, m`
     - Bus indices
   * - :math:`g`
     - Generator index
   * - :math:`l`
     - Transmission line index

Parameters
~~~~~~~~~~
.. list-table::
   :widths: 20 80
   :header-rows: 0

   * - :math:`S_{\text{base}}`
     - Base power [MVA]
   * - :math:`P^D_n`
     - Active power demand at bus :math:`n` [MW]
   * - :math:`\theta_{\text{max}}, \theta_{\text{min}}`
     - Voltage angle limits [radians]
   * - :math:`\overline{P}_g, \underline{P}_g`
     - Generation limits for :math:`g` [MW]
   * - :math:`c_{2,g}, c_{1,g}, c_{0,g}`
     - Cost coefficients for :math:`g`
   * - :math:`x_l`
     - Reactance of line :math:`l` [p.u.]
   * - :math:`\overline{P}_l`
     - Line capacity for :math:`l` [MW]
   * - :math:`C^F_{l,n}, C^T_{l,n}`
     - Incidence matrices for "from" and "to" buses of line :math:`l`
   * - :math:`C^G_{g,n}`
     - Incidence matrix for generator :math:`g` at bus :math:`n`

Variables
~~~~~~~~~
.. list-table::
   :widths: 20 80
   :header-rows: 0

   * - :math:`P_g`
     - Active power generation of :math:`g` [MW]
   * - :math:`\theta_n`
     - Voltage angle at bus :math:`n` [radians]
   * - :math:`P^F_l, P^T_l`
     - Power flow on line :math:`l` (from/to) [MW]

DC Optimal Power Flow (DC-OPF)
------------------------------

The DC-OPF formulation is a simplified version of the AC-OPF that assumes:

- Voltage magnitudes are fixed at 1.0 p.u.
- Reactive power flows are ignored.
- Small angle differences between buses.

Objective Function
~~~~~~~~~~~~~~~~~
The objective is to minimize the total generation cost:

.. math::
   \text{Minimize: } \sum_{g} \left( c_{2,g} \cdot P_g^2 + c_{1,g} \cdot P_g + c_{0,g} \right)

Constraints
~~~~~~~~~~~

1. **Active Power Balance**:

   .. math::
      \sum_{g} C^G_{g,n} \cdot P_g - P^D_n = \sum_{l} \left( C^F_{l,n} \cdot P^F_l + C^T_{l,n} \cdot P^T_l \right), \quad \forall n

2. **Line Flow Equations**:

   .. math::
      P^F_l = \frac{1}{x_l} \cdot \sum_{n} \left( C^F_{l,n} - C^T_{l,n} \right) \cdot \theta_n, \quad \forall l

   .. math::
      P^T_l = \frac{1}{x_l} \cdot \sum_{n} \left( C^T_{l,n} - C^F_{l,n} \right) \cdot \theta_n, \quad \forall l

3. **Line Flow Limits**:

   .. math::
      -\overline{P}_l \leq P^F_l, P^T_l \leq \overline{P}_l, \quad \forall l

4. **Generator Limits**:

   .. math::
      \underline{P}_g \leq P_g \leq \overline{P}_g, \quad \forall g

5. **Voltage Angle Reference**:

   .. math::
      \theta_{\text{slack}} = 0
