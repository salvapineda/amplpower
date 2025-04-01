------------------------------
AC Optimal Power Flow (rectangular coordinates)
------------------------------

This formulation uses rectangular coordinates for voltage representation, where the voltage at node :math:`n` is represented as:

.. math::
   \mathcal{V}_n = e_n + j f_n

where :math:`e_n` and :math:`f_n` are the real and imaginary parts of the voltage, respectively. The AC-OPF formulation considers the full nonlinear power flow equations, including both active and reactive power, and the real and imaginary parts of the voltage.
The AC-OPF problem is formulated as a non-convex optimization problem. In particular, the AC_OPF is formulated as a non-convex quadratically constrainted quadratic program (QCQP).

Objective Function
------------------
The objective is to minimize the total generation cost:

.. math::
   \text{Minimize: } \sum_{g} \left( a_{g} \cdot p_g^2 + b_{g} \cdot p_g + c_{g} \right)

Constraints
-----------

1. **Power Balance Equations**:

   .. math::
      \sum_{g\in\mathcal{G}_n} p_g - p^d_n = G^{sh}_n v_n^{(2)} + \sum_{l} \left( F_{ln} \cdot p^f_l + T_{ln} \cdot p^t_l \right), \quad \forall n

   .. math::
      \sum_{g\in\mathcal{G}_n} q_g - q^d_n = -B^{sh}_n v_n^{(2)} + \sum_{l} \left( F_{ln} \cdot q^f_l + T_{ln} \cdot q^t_l \right), \quad \forall n

2. **Line Flow Equations**:

   .. math::
      p^f_l = G^{ff}_l v_n^{(2)} + G^{ft}_l c^{ft}_l + B^{ft}_l s^{ft}_l, \quad \forall (l, n, m): F_{ln} = 1, T_{lm} = 1

   .. math::
      q^f_l = -B^{ff}_l v_n^{(2)} + G^{ft}_l s^{ft}_l - B^{ft}_l c^{ft}_l, \quad \forall (l, n, m): F_{ln} = 1, T_{lm} = 1

   .. math::
      p^t_l = G^{tt}_l v_m^{(2)} + G^{tf}_l c^{ft}_l + B^{tf}_l s^{ft}_l, \quad \forall (l, n, m): F_{ln} = 1, T_{lm} = 1

   .. math::
      q^t_l = -B^{tt}_l v_m^{(2)} + G^{tf}_l s^{ft}_l - B^{tf}_l c^{ft}_l, \quad \forall (l, n, m): F_{ln} = 1, T_{lm} = 1

3. **Rectangular Definitions**:

   .. math::
      v_{n}^{(2)} = e_n^2 + f_n^2, \quad \forall n

   .. math::
      c^{ft}_l = e_n e_m + f_n f_m, \quad \forall (l,n,m): F_{ln} = 1, T_{lm} = 1

   .. math::
      s^{ft}_l = f_n e_m - e_n f_m, \quad \forall (l,n,m): F_{ln} = 1, T_{lm} = 1

4. **Generator Limits**:

   .. math::
      \underline{p}_g \leq p_g \leq \overline{p}_g, \quad \forall g

   .. math::
      \underline{q}_g \leq q_g \leq \overline{q}_g, \quad \forall g

5. **Line Flow Limits**:

   .. math::
      (p^f_l)^2 + (q^f_l)^2 \leq (\overline{s}_l)^2, \quad \forall l

   .. math::
      (p^t_l)^2 + (q^t_l)^2 \leq (\overline{s}_l)^2, \quad \forall l

6. **Voltage Magnitude Limits**:

   .. math::
      (\underline{v}_n)^2 \leq v_n^{(2)} \leq (\overline{v}_n)^2, \quad \forall n

7. **Voltage Angle Reference**:

   .. math::
      f_0 = 0
