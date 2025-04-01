------------------------------
AC Optimal Power Flow (polar coordinates)
------------------------------

This formulation uses polar coordinates for voltage representation, where the voltage at node :math:`n` is represented as:

.. math::
   \mathcal{V}_n = v_n e^{j\theta_n}

where :math:`v_n` is the voltage magnitude and :math:`\theta_n` is the voltage angle. The AC-OPF formulation considers the full nonlinear power flow equations, including both active and reactive power, voltage magnitudes, and angles. Mathematically, the AC-OPF problem is formulated as a non-convex optimization problem.

Objective Function
------------------
The objective is to minimize the total generation cost:

.. math::
   \text{Minimize: } \sum_{g} \left( a_{g} \cdot p_g^2 + b_{g} \cdot p_g + c_{g} \right)

Constraints
-----------

1. **Power Balance Equations**:

   .. math::
      \sum_{g\in\mathcal{G}_n} p_g - p^d_n = v_n^2 G^{sh}_n + \sum_{l} \left( F_{ln} \cdot p^f_l + T_{ln} \cdot p^t_l \right), \quad \forall n

   .. math::
      \sum_{g\in\mathcal{G}_n} q_g - q^d_n = -v_n^2 G^{sh}_n + \sum_{l} \left( F_{ln} \cdot q^f_l + T_{ln} \cdot q^t_l \right), \quad \forall n

2. **Line Flow Equations**:

   .. math::
      p^f_l = v_n^2 G^{ff}_l + v_n v_m \left( G^{ft}_l \cos(\theta_{nm}) + B^{ft}_l \sin(\theta_{nm}) \right), \quad \forall (l, n, m): F_{ln} = 1, T_{lm} = 1

   .. math::
      q^f_l = -v_n^2 B^{ff}_l + v_n v_m \left( G^{ft}_l \sin(\theta_{nm}) - B^{ft}_l \cos(\theta_{nm}) \right), \quad \forall (l, n, m): F_{ln} = 1, T_{lm} = 1

   .. math::
      p^t_l = v_m^2 G^{tt}_l + v_n v_m \left( G^{tf}_l \cos(\theta_{mn}) + B^{tf}_l \sin(\theta_{mn}) \right), \quad \forall (l, n, m): F_{ln} = 1, T_{lm} = 1

   .. math::
      q^t_l = -v_m^2 B^{tt}_l + v_n v_m \left( G^{tf}_l \sin(\theta_{mn}) - B^{tf}_l \cos(\theta_{mn}) \right), \quad \forall (l, n, m): F_{ln} = 1, T_{lm} = 1

3. **Generator Limits**:

   .. math::
      \underline{p}_g \leq p_g \leq \overline{p}_g, \quad \forall g

   .. math::
      \underline{q}_g \leq q_g \leq \overline{q}_g, \quad \forall g

4. **Line Flow Limits**:

   .. math::
      (p^f_l)^2 + (q^f_l)^2 \leq (\overline{s}_l)^2, \quad \forall l

   .. math::
      (p^t_l)^2 + (q^t_l)^2 \leq (\overline{s}_l)^2, \quad \forall l

5. **Voltage Magnitude Limits**:

   .. math::
      \underline{v}_n \leq v_n \leq \overline{v}_n, \quad \forall n

   .. math::
      \underline{\theta}_n \leq \theta_n \leq \overline{\theta}_n, \quad \forall n

6. **Voltage Angle Reference**:

   .. math::
      \theta_{0} = 0
