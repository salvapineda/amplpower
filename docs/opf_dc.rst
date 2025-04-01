------------------------------
DC Optimal Power Flow
------------------------------

The DC-OPF formulation is a simplified version of the AC-OPF that assumes:

- Voltage magnitudes are fixed at 1.0 p.u.
- Reactive power flows are ignored.
- Small angle differences between buses.

If the generating cost functions for all generators are linear, the DC-OPF formulation becomes a linear programming (LP) problem. Otherwise, it is a convex quadratic programming (QP) problem due to the quadratic term in the objective function.

Objective Function
------------------
The objective is to minimize the total generation cost:

.. math::
   \text{Minimize: } \sum_{g} \left( a_g \cdot p_g^2 + b_{g} \cdot p_g + c_{g} \right)

Constraints
-----------

1. **Active Power Balance**:

   .. math::
      \sum_{g\in\mathcal{G}_n} p_g - p^d_n = \sum_{l} \left( F_{ln} \cdot p^f_l + T_{ln} \cdot p^t_l \right), \quad \forall n

2. **Line Flow Equations**:

   .. math::
      p^f_l = \frac{1}{x_l} \cdot \left( \theta_n - \theta_m\right), \quad \forall (l,n,m): F_{l,n} = 1, T_{l,m} = 1

   .. math::
      p^t_l = \frac{1}{x_l} \cdot \left( \theta_m - \theta_n \right), \quad \forall (l,n,m): F_{l,n} = 1, T_{l,m} = 1

3. **Line Flow Limits**:

   .. math::
      -\overline{s}_l \leq p^f_l \leq \overline{s}_l, \quad \forall l
   .. math::
      -\overline{s}_l \leq p^t_l \leq \overline{s}_l, \quad \forall l

4. **Generator Limits**:

   .. math::
      \underline{p}_g \leq p_g \leq \overline{p}_g, \quad \forall g

5. **Voltage Magnitude Limits**:

   .. math::
      \underline{\theta}_n \leq \theta_n \leq \overline{\theta}_n, \quad \forall n

6. **Voltage Angle Reference**:

   .. math::
      \theta_{0} = 0
