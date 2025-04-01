===========================
Nomenclature
===========================

This section provides an overview of the formulations and optimization models used in the package.

Indexes
-------
.. list-table::
   :widths: 20 80
   :header-rows: 0

   * - :math:`g`
     - Generator index
   * - :math:`\mathcal{G}_n`
     - Set of generators at bus :math:`n`
   * - :math:`n, m`
     - Bus indices
   * - :math:`l`
     - Transmission line index

Parameters
----------
.. list-table::
   :widths: 20 80
   :header-rows: 0

   * - :math:`S_{\text{base}}`
     - Base power [MVA]
   * - :math:`a_{g}, b_{g}, c_{g}`
     - Cost coefficients for :math:`g`
   * - :math:`\overline{p}_g, \underline{p}_g`
     - Active power generation limits for :math:`g` [MW]
   * - :math:`\overline{q}_g, \underline{q}_g`
     - Reactive power generation limits for :math:`g` [MW]
   * - :math:`p^d_n, q^d_n`
     - Active and reactive power demand at node :math:`n` [MW, MVAr]
   * - :math:`\overline{\theta}_n, \underline{\theta}_n`
     - Voltage angle limits at bus :math:`n` [radians]
   * - :math:`\underline{v}_n, \overline{v}_n`
     - Voltage magnitude limits at node :math:`n` [p.u.]
   * - :math:`G^{sh}_n, B^{sh}_n`
     - Shunt conductance and susceptance at node :math:`n` [p.u.]
   * - :math:`r_l, x_l`
     - Series resistance and reactance of branch :math:`l` [p.u.]
   * - :math:`b_l`
     - Total charging susceptance of branch :math:`l` [p.u.]
   * - :math:`\tau_l, \gamma_l`
     - Tap ratio magnitude and phase shift angle of branch :math:`l`
   * - :math:`\overline{s}_l`
     - Apparent power flow limit for branch :math:`l` [MVA]
   * - :math:`F_{l,n}, T_{l,n}`
     - Incidence matrices for "from" and "to" buses of line :math:`l`

We also define the following parameters for the branch :math:`l`:

.. math::
   \begin{align}
   & G^{ff}_l + jB^{ff}_l = \left(\tfrac{1}{r_l+jx_l}+j\tfrac{b_l}{2} \right)\tfrac{1}{\tau_l^2} \\
   & G^{ft}_l + jB^{ft}_l = -\tfrac{1}{r_l+jx_l}\tfrac{1}{\tau_l e^{-j\gamma_l}} \\
   & G^{tf}_l + jB^{tf}_l = -\tfrac{1}{r_l+jx_l}\tfrac{1}{\tau_l e^{j\gamma_l}} \\
   & G^{tt}_l + jB^{tt}_l = \left(\tfrac{1}{r_l+jx_l}+j\tfrac{b_l}{2} \right)
   \end{align}

Variables
---------
.. list-table::
   :widths: 20 80
   :header-rows: 0

   * - :math:`p_g`
     - Active power generation of :math:`g` [MW]
   * - :math:`q_g`
     - Reactive power generation of :math:`g` [MVAr]
   * - :math:`v_n`
     - Voltage magnitude at node :math:`n` [p.u.]
   * - :math:`\theta_n`
     - Voltage angle at bus :math:`n` [radians]
   * - :math:`v_n^{(2)}`
     -  Voltage magnitude squared at node :math:`n` [p.u.]
   * - :math:`e_n, f_n`
     - Real and imaginary parts of voltage at node :math:`n` [p.u.]
   * - :math:`p^f_l, p^t_l`
     - Active power flow on branch :math:`l` (from/to) [MW]
   * - :math:`q^f_l, q^t_l`
     - Reactive power flow on branch :math:`l` (from/to) [MVAr]
