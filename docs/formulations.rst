===========================
Formulations and Models
===========================

This section provides an overview of the formulations and optimization models used in the package.

Nomenclature
------------

### Indexes
- **\( i, j \)**: Buses.
- **\( g \)**: Generators.
- **\( l \)**: Transmission lines.

### Parameters
- **\( S_{\text{base}} \)**: Base power in MVA.
- **\( P_{D,i} \)**: Active power demand at bus \( i \) (MW).
- **\( \theta_{\text{max}}, \theta_{\text{min}} \)**: Voltage angle limits (radians).
- **\( P_{g,\text{max}}, P_{g,\text{min}} \)**: Active power generation limits for generator \( g \) (MW).
- **\( c_{2,g}, c_{1,g}, c_{0,g} \)**: Quadratic, linear, and constant cost coefficients for generator \( g \).
- **\( x_l \)**: Reactance of line \( l \) (p.u.).
- **\( f_l, t_l \)**: From and to buses of line \( l \).

### Variables
- **\( P_g \)**: Active power generation of generator \( g \) (MW).
- **\( \theta_i \)**: Voltage angle at bus \( i \) (radians).
- **\( P_{f,l}, P_{t,l} \)**: Active power flow on line \( l \) from the "from" and "to" ends (MW).

Formulations
------------

- **P**: Active power (MW)
- **Q**: Reactive power (MVAR)
- **V**: Voltage magnitude (p.u.)
- **θ**: Voltage angle (radians)
- **Y**: Admittance matrix

Optimization Models
-------------------

1. **AC Optimal Power Flow (AC-OPF)**:
   - Solves the full nonlinear power flow equations.
   - Objective: Minimize generation cost or losses.

2. **DC Optimal Power Flow (DC-OPF)**:
   - Linearized version of the AC-OPF.
   - Assumes small angle differences and constant voltage magnitudes.

DC Optimal Power Flow (DC-OPF)
------------------------------

The DC-OPF formulation is a simplified version of the AC-OPF that assumes:
- Voltage magnitudes are fixed at 1.0 p.u.
- Reactive power flows are ignored.
- Small angle differences between buses.

### Objective Function
The objective is to minimize the total generation cost:
\[
\text{Minimize: } \sum_{g} \left( c_{2,g} \cdot P_g^2 + c_{1,g} \cdot P_g + c_{0,g} \right)
\]

### Constraints

1. **Active Power Balance**:
   For each bus \( i \):
   \[
   \sum_{g} P_g \cdot \delta_{g,i} - P_{D,i} = \sum_{l} \left( P_{f,l} \cdot \delta_{f_l,i} + P_{t,l} \cdot \delta_{t_l,i} \right)
   \]
   where \( \delta_{g,i} \), \( \delta_{f_l,i} \), and \( \delta_{t_l,i} \) are incidence indicators.

2. **Line Flow Equations**:
   For each line \( l \):
   \[
   P_{f,l} = \frac{1}{x_l} \cdot \left( \theta_{f_l} - \theta_{t_l} \right)
   \]
   \[
   P_{t,l} = \frac{1}{x_l} \cdot \left( \theta_{t_l} - \theta_{f_l} \right)
   \]

3. **Generator Limits**:
   For each generator \( g \):
   \[
   P_{g,\text{min}} \leq P_g \leq P_{g,\text{max}}
   \]

4. **Voltage Angle Reference**:
   The voltage angle at the slack bus is fixed:
   \[
   \theta_{\text{slack}} = 0
   \]

3. **Security-Constrained OPF (SC-OPF)**:
   - Includes contingency constraints for system reliability.
   - Ensures the system operates securely under predefined contingencies.

Feel free to expand this section with additional formulations or constraints as needed.
