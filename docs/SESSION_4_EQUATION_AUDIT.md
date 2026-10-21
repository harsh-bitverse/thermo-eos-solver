# CHL2002 Thermodynamics Project — Session 4 Equation Audit
## Residual Properties, Fugacity, Reference States, and Phase Equilibrium

**Document Version:** 1.0  
**Date:** October 21, 2026  
**Status:** Approved Pre-Implementation Mathematical Specification  

---

## 1. Governing Thermodynamic Equations and Units

All calculations in the core engine use strict SI units:
* Pressure $P$: $[\text{Pa}]$
* Temperature $T$: $[\text{K}]$
* Molar volume $v$: $[\text{m}^3/\text{mol}]$
* Universal gas constant $R = 8.31446261815324 \text{ J}/(\text{mol}\cdot\text{K})$
* PR attractive parameter $a(T)$: $[\text{Pa}\cdot\text{m}^6/\text{mol}^2]$
* Co-volume parameter $b$: $[\text{m}^3/\text{mol}]$
* Temperature derivative $\frac{da}{dT}$: $[\text{Pa}\cdot\text{m}^6/(\text{mol}^2\cdot\text{K})]$
* Enthalpy ($H^{\text{res}}, H$): $[\text{J}/\text{mol}]$
* Entropy ($S^{\text{res}}, S$): $[\text{J}/(\text{mol}\cdot\text{K})]$
* Fugacity coefficient $\phi$: dimensionless $[-]$
* Fugacity $f = \phi P$: $[\text{Pa}]$

### 1.1 Logarithmic Departure Term $L(Z, B)$
For any physically admissible root ($Z > B \iff v > b$):
$$L = \ln \left( \frac{Z + (1 + \sqrt{2})B}{Z + (1 - \sqrt{2})B} \right)$$
where $B = \frac{bP}{RT}$ and $A = \frac{a(T)P}{(RT)^2}$.

### 1.2 Residual Enthalpy $H^{\text{res}}$
$$H^{\text{res}} = H(T, P) - H^{\text{ig}}(T, P) = RT(Z - 1) + \frac{T \frac{da}{dT} - a(T)}{2\sqrt{2}b} L \quad [\text{J}/\text{mol}]$$

### 1.3 Residual Entropy $S^{\text{res}}$
$$S^{\text{res}} = S(T, P) - S^{\text{ig}}(T, P) = R \ln(Z - B) + \frac{\frac{da}{dT}}{2\sqrt{2}b} L \quad [\text{J}/(\text{mol}\cdot\text{K})]$$

### 1.4 Pure-Component Fugacity Coefficient $\phi$ and Fugacity $f$
$$\ln \phi = (Z - 1) - \ln(Z - B) - \frac{A}{2\sqrt{2}B} L$$
$$\phi = \exp(\ln \phi), \qquad f = \phi P \quad [\text{Pa}]$$

### 1.5 Rigorous Entropy Identity
Combining $H^{\text{res}}$, $S^{\text{res}}$, and $\ln \phi$:
$$\frac{H^{\text{res}}}{T} - R \ln \phi = R(Z - 1) + \frac{T \frac{da}{dT} - a}{2\sqrt{2}b T} L - R(Z - 1) + R\ln(Z - B) + \frac{RA}{2\sqrt{2}B} L$$
Since $\frac{RA}{B} = \frac{a(T)}{b T}$, the attractive parameter $a(T)$ cancels identically:
$$= R \ln(Z - B) + \frac{\frac{da}{dT}}{2\sqrt{2}b} L \equiv S^{\text{res}}$$
Hence, the Maxwell identity:
$$S^{\text{res}} = \frac{H^{\text{res}}}{T} - R \ln \phi$$
holds with analytical exactness and machine precision ($\le 10^{-15}\text{ J}/(\text{mol}\cdot\text{K})$).

---

## 2. Low-Pressure Numerical Cancellation Control

At low pressure ($P \to 0$), $B \to 0$ and $Z \to 1$. Evaluating $L / (2\sqrt{2}B)$ and $(Z - 1) - \ln(Z - B)$ directly can induce severe loss of precision through floating-point cancellation.

### 2.1 Regularization of $L / (2\sqrt{2}B)$
Let:
$$u = \frac{2\sqrt{2}B}{Z + (1 - \sqrt{2})B}$$
Then the argument of the logarithm satisfies:
$$\frac{Z + (1 + \sqrt{2})B}{Z + (1 - \sqrt{2})B} = 1 + u$$
Thus $L = \ln(1 + u)$.
Define the dimensionless ratio $I_L = \frac{L}{2\sqrt{2}B}$:
$$I_L = \frac{\ln(1 + u)}{2\sqrt{2}B} = \frac{1}{Z + (1 - \sqrt{2})B} \cdot \frac{\ln(1 + u)}{u}$$

* **For $|u| \ge 10^{-5}$:** Evaluate via `std::log1p`:
  $$L = \text{std::log1p}(u), \qquad I_L = \frac{L}{2\sqrt{2}B}$$
* **For $|u| < 10^{-5}$:** Use the 4-term Taylor expansion of $\frac{\ln(1+u)}{u}$:
  $$\frac{\ln(1 + u)}{u} = 1 - \frac{u}{2} + \frac{u^2}{3} - \frac{u^3}{4} + \mathcal{O}(u^4)$$
  $$I_L = \frac{1}{Z + (1 - \sqrt{2})B} \left( 1 - \frac{u}{2} + \frac{u^2}{3} - \frac{u^3}{4} \right), \qquad L = 2\sqrt{2}B \cdot I_L$$
* **Tracer continuity:** At $|u| = 10^{-5}$, truncation error is $\frac{u^4}{5} = 2 \times 10^{-21} \ll \epsilon_{\text{mach}} \approx 2.22 \times 10^{-16}$. Direct evaluation and Taylor series match to within floating-point roundoff.

Furthermore, in $H^{\text{res}}$ and $S^{\text{res}}$:
$$\frac{L}{2\sqrt{2}b} = I_L \cdot \frac{P}{RT}$$
which smoothly approaches zero as $P \to 0$ without any division by $b$ cancellation.

### 2.2 Regularization of $(Z - 1) - \ln(Z - B)$
Let $\epsilon = Z - 1 - B$. Then $Z - B = 1 + \epsilon$.
Using $\ln(1 + \epsilon) = \epsilon - \frac{\epsilon^2}{2} + \frac{\epsilon^3}{3} - \frac{\epsilon^4}{4} + \dots$:
$$(Z - 1) - \ln(Z - B) = (Z - 1) - \left( \epsilon - \frac{\epsilon^2}{2} + \frac{\epsilon^3}{3} - \frac{\epsilon^4}{4} \right) = B + \epsilon^2 \left( \frac{1}{2} - \frac{\epsilon}{3} + \frac{\epsilon^2}{4} \right)$$
* **For $|\epsilon| \ge 10^{-5}$:** $(Z - 1.0) - \text{std::log}(Z - B)$.
* **For $|\epsilon| < 10^{-5}$:** $B + \epsilon^2 \left( \frac{1}{2} - \frac{\epsilon}{3} + \frac{\epsilon^2}{4} \right)$.
This eliminates catastrophic cancellation between $Z - 1$ and $\ln(Z - B)$ down to $P = 0$.

---

## 3. Reference-State Convention for Total $H$ and $S$

### 3.1 Standard Reference State
The reference state is defined as:
$$T_{\mathrm{ref}} = 298.15\text{ K} \quad (25.0^\circ\text{C}), \qquad P_{\mathrm{ref}} = 100\,000.0\text{ Pa} \quad (1.0\text{ bar})$$
At this state, the reference convention dictates:
$$H(T_{\mathrm{ref}}, P_{\mathrm{ref}}) = 0, \qquad S(T_{\mathrm{ref}}, P_{\mathrm{ref}}) = 0$$

### 3.2 Ideal-Gas Integrals via NIST Shomate
For $t = T / 1000$ and $t_0 = T_{\mathrm{ref}} / 1000$:
$$\Delta H^{\text{ig}}(T_{\mathrm{ref}} \to T) = \int_{T_{\mathrm{ref}}}^{T} C_p^\circ(T') dT' = 1000 \left[ A(t - t_0) + \frac{B}{2}(t^2 - t_0^2) + \frac{C}{3}(t^3 - t_0^3) + \frac{D}{4}(t^4 - t_0^4) - E\left(\frac{1}{t} - \frac{1}{t_0}\right) \right]$$
$$\Delta S^{\text{ig}}(T_{\mathrm{ref}} \to T) = \int_{T_{\mathrm{ref}}}^{T} \frac{C_p^\circ(T')}{T'} dT' = A \ln\left(\frac{t}{t_0}\right) + B(t - t_0) + \frac{C}{2}(t^2 - t_0^2) + \frac{D}{3}(t^3 - t_0^3) - \frac{E}{2}\left(\frac{1}{t^2} - \frac{1}{t_0^2}\right)$$

### 3.3 Reference Departure Corrections
To guarantee that $H(T_{\mathrm{ref}}, P_{\mathrm{ref}}) = 0$ and $S(T_{\mathrm{ref}}, P_{\mathrm{ref}}) = 0$:
$$H(T, P) = \Delta H^{\text{ig}}(T_{\mathrm{ref}} \to T) + H^{\text{res}}(T, P) - H^{\text{res}}(T_{\mathrm{ref}}, P_{\mathrm{ref}})$$
$$S(T, P) = \Delta S^{\text{ig}}(T_{\mathrm{ref}} \to T) - R \ln\left(\frac{P}{P_{\mathrm{ref}}}\right) + S^{\text{res}}(T, P) - S^{\text{res}}(T_{\mathrm{ref}}, P_{\mathrm{ref}})$$

* **Reference Root Selection:** At $(T_{\mathrm{ref}}, P_{\mathrm{ref}})$, the reference departures $H^{\text{res}}(T_{\mathrm{ref}}, P_{\mathrm{ref}})$ and $S^{\text{res}}(T_{\mathrm{ref}}, P_{\mathrm{ref}})$ are evaluated on the **vapor-like root** (`vapor_root()`), corresponding to the physical ideal-gas asymptotic fluid branch.
* **Calorimetric Distinction:** These total properties represent state increments relative to the $(T_{\mathrm{ref}}, P_{\mathrm{ref}})$ standard state and must not be confused with absolute enthalpy of formation or third-law calorimetric entropy.

---

## 4. Multiple Roots and Thermodynamic Phase Preference

When multiple physical roots exist at $(T, P)$ ($Z_{\min} < Z_{\text{mid}} < Z_{\max}$):
1. **Mechanical Stability Filter:**
   * The intermediate root $Z_{\text{mid}}$ has $\left(\frac{\partial P}{\partial v}\right)_T > 0$. It is mechanically unstable and strictly excluded from phase equilibrium comparisons.
2. **Gibbs Free Energy Comparison:**
   For a pure component, the chemical potentials of liquid and vapor satisfy:
   $$G_{\text{liq}} - G_{\text{vap}} = G^{\text{res}}_{\text{liq}} - G^{\text{res}}_{\text{vap}} = RT \ln \left( \frac{\phi_{\text{liq}}}{\phi_{\text{vap}}} \right)$$
   Define $\Delta G_{\text{res}} = G^{\text{res}}_{\text{liq}} - G^{\text{res}}_{\text{vap}}$.
3. **Classification Criteria:**
   * If $\ln(\phi_{\text{liq}} / \phi_{\text{vap}}) < -\text{tol}$: $\Delta G_{\text{res}} < 0 \implies G_{\text{liq}} < G_{\text{vap}}$. **Liquid is thermodynamically preferred** (subcooled liquid).
   * If $\ln(\phi_{\text{liq}} / \phi_{\text{vap}}) > +\text{tol}$: $\Delta G_{\text{res}} > 0 \implies G_{\text{vap}} < G_{\text{liq}}$. **Vapor is thermodynamically preferred** (superheated vapor).
   * If $|\ln(\phi_{\text{liq}} / \phi_{\text{vap}})| \le \text{tol}$: **Two-phase equilibrium (saturation state)** within numerical tolerance.
4. **Tolerance:**
   We enforce `EQUILIBRIUM_FUGACITY_TOL = 1.0e-5` on $|\phi_{\text{liq}} - \phi_{\text{vap}}| / \phi_{\text{vap}}$.
5. **Single Physical Root:**
   When only one physical root is found, the system reports `PhasePreference::SINGLE_PHASE` without claiming phase coexistence.

---

## 5. Component Database Limitations & Temperature Validity Ranges

| Component | Formula | Valid Range $T$ (K) | $T_{\mathrm{ref}} = 298.15\text{ K}$ Status | Handling / Policy |
| :--- | :--- | :--- | :--- | :--- |
| **Methane** | $\text{CH}_4$ | $298.0 - 1300.0$ | Valid ($\in [298.0, 1300.0]$) | Full $H, S$ support. |
| **Ethane** | $\text{C}_2\text{H}_6$ | $298.0 - 1200.0$ | Valid ($\in [298.0, 1200.0]$) | Full $H, S$ support. |
| **Propane** | $\text{C}_3\text{H}_8$ | $298.0 - 1200.0$ | Valid ($\in [298.0, 1200.0]$) | Full $H, S$ support. |
| **n-Butane** | $\text{C}_4\text{H}_{10}$ | $298.0 - 1200.0$ | Valid ($\in [298.0, 1200.0]$) | Full $H, S$ support. |
| **Nitrogen** | $\text{N}_2$ | $100.0 - 500.0$ | Valid ($\in [100.0, 500.0]$) | Full $H, S$ support. |
| **Oxygen** | $\text{O}_2$ | $100.0 - 700.0$ | Valid ($\in [100.0, 700.0]$) | Full $H, S$ support. |
| **Carbon Dioxide** | $\text{CO}_2$ | $298.0 - 1200.0$ | Valid ($\in [298.0, 1200.0]$) | Full $H, S$ support. |
| **Water Vapor** | $\text{H}_2\text{O}$ | $500.0 - 1700.0$ | **Invalid** ($T_{\mathrm{ref}} < 500.0\text{ K}$) | Gas-phase Shomate parameters begin at $500\text{ K}$. Total $H(T, P)$ and $S(T, P)$ referencing $298.15\text{ K}$ will explicitly throw `std::out_of_range`. Residual properties ($H^{\text{res}}, S^{\text{res}}, \phi, f$) remain fully available. |

---

## 6. Independently Calculated Benchmark Targets

The following benchmark values are independently derived using the PR-1976 model and authenticated NIST Shomate parameters:

### 6.1 Methane ($\text{CH}_4$) at $300\text{ K}, 5\text{ MPa}$
* $Z = 0.901844$, $v = 4.4990 \times 10^{-4}\text{ m}^3/\text{mol}$
* $H^{\text{res}} = -902.3382\text{ J/mol}$
* $S^{\text{res}} = -2.144565\text{ J}/(\text{mol}\cdot\text{K})$
* $\ln \phi = -0.103823 \implies \phi = 0.901385$, $f = 4.506926\text{ MPa}$
* Identity error: $S^{\text{res}} - (H^{\text{res}}/T - R \ln \phi) = 0.00 \times 10^{-16}\text{ J}/(\text{mol}\cdot\text{K})$
* Reference state ($298.15\text{ K}, 1\text{ bar}$): $Z = 0.997783$, $H^{\text{res}}_{\text{ref}} = -18.1480\text{ J/mol}$, $S^{\text{res}}_{\text{ref}} = -0.042423\text{ J}/(\text{mol}\cdot\text{K})$
* $\Delta H^{\text{ig}} = 66.0074\text{ J/mol}$, $\Delta S^{\text{ig}} = 0.220706\text{ J}/(\text{mol}\cdot\text{K})$
* Total $H = -818.1829\text{ J/mol}$, Total $S = -34.407805\text{ J}/(\text{mol}\cdot\text{K})$

### 6.2 Carbon Dioxide ($\text{CO}_2$) at $300\text{ K}, 4\text{ MPa}$
* $Z = 0.752334$, $v = 4.6914 \times 10^{-4}\text{ m}^3/\text{mol}$
* $H^{\text{res}} = -1990.9877\text{ J/mol}$
* $S^{\text{res}} = -4.734694\text{ J}/(\text{mol}\cdot\text{K})$
* $\ln \phi = -0.228750 \implies \phi = 0.795528$, $f = 3.182110\text{ MPa}$
* Identity error: $8.88 \times 10^{-16}\text{ J}/(\text{mol}\cdot\text{K})$
* Reference state ($298.15\text{ K}, 1\text{ bar}$): $Z = 0.994521$, $H^{\text{res}}_{\text{ref}} = -41.1744\text{ J/mol}$, $S^{\text{res}}_{\text{ref}} = -0.092609\text{ J}/(\text{mol}\cdot\text{K})$
* $\Delta H^{\text{ig}} = 68.7716\text{ J/mol}$, $\Delta S^{\text{ig}} = 0.229948\text{ J}/(\text{mol}\cdot\text{K})$
* Total $H = -1881.0417\text{ J/mol}$, Total $S = -35.083187\text{ J}/(\text{mol}\cdot\text{K})$

### 6.3 Carbon Dioxide ($\text{CO}_2$) at $300\text{ K}, 7\text{ MPa}$
* $Z = 0.198057$, $v = 7.0575 \times 10^{-5}\text{ m}^3/\text{mol}$ (Single liquid root)
* $H^{\text{res}} = -9374.3721\text{ J/mol}$
* $S^{\text{res}} = -27.514960\text{ J}/(\text{mol}\cdot\text{K})$
* $\ln \phi = -0.448970 \implies \phi = 0.638285$, $f = 4.467995\text{ MPa}$
* Total $H = -9264.4261\text{ J/mol}$, Total $S = -62.516357\text{ J}/(\text{mol}\cdot\text{K})$

### 6.4 Propane ($\text{C}_3\text{H}_8$) at $320\text{ K}, 1\text{ MPa}$ (Three Roots)
* Reference state ($298.15\text{ K}, 1\text{ bar}$, vapor root): $Z = 0.983423$, $H^{\text{res}}_{\text{ref}} = -111.9835\text{ J/mol}$, $S^{\text{res}}_{\text{ref}} = -0.238539\text{ J}/(\text{mol}\cdot\text{K})$
* $\Delta H^{\text{ig}} = 1848.1478\text{ J/mol}$, $\Delta S^{\text{ig}} = 5.980242\text{ J}/(\text{mol}\cdot\text{K})$
* **Root 0 (Liquid):** $Z = 0.036453$, $H^{\text{res}} = -14854.34\text{ J/mol}$, $S^{\text{res}} = -48.2035\text{ J}/(\text{mol}\cdot\text{K})$, $\phi = 1.239280$, Total $H = -12894.21\text{ J/mol}$, Total $S = -61.13\text{ J}/(\text{mol}\cdot\text{K})$
* **Root 1 (Intermediate):** $Z = 0.091826$, $H^{\text{res}} = -8210.60\text{ J/mol}$, $S^{\text{res}} = -28.3698\text{ J}/(\text{mol}\cdot\text{K})$, $\phi = 1.385615$, UNSTABLE
* **Root 2 (Vapor):** $Z = 0.850555$, $H^{\text{res}} = -1131.04\text{ J/mol}$, $S^{\text{res}} = -2.3575\text{ J}/(\text{mol}\cdot\text{K})$, $\phi = 0.868006$, Total $H = +829.09\text{ J/mol}$, Total $S = -15.28\text{ J}/(\text{mol}\cdot\text{K})$
* **Phase Preference:** $\phi_{\text{vap}} = 0.868006 < \phi_{\text{liq}} = 1.239280 \implies \ln(\phi_{\text{liq}}/\phi_{\text{vap}}) = 0.3559 > 0 \implies$ **Vapor phase preferred** (superheated).

### 6.5 Propane Saturation Benchmark ($320\text{ K}$)
* $P^{\text{sat}} = 1\,603\,614.23\text{ Pa}$ ($1.603614\text{ MPa}$)
* $\phi_{\text{liq}} = 0.78988397$, $\phi_{\text{vap}} = 0.78988397$
* Relative difference: $\frac{|\phi_{\text{liq}} - \phi_{\text{vap}}|}{\phi_{\text{vap}}} \approx 1.09 \times 10^{-13}$
* $Z_{\text{liq}} = 0.057709$, $Z_{\text{vap}} = 0.736416$
* Result: `PhasePreference::EQUILIBRIUM`
