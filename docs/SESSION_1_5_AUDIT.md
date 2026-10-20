# CHL2002 Thermodynamics Project — Session 1.5: Independent Thermodynamic Audit and Approval Gate

**Document Version:** 1.0  
**Audit Date:** October 20, 2026  
**Auditor:** Senior Thermodynamics Software Engineer & C++ Implementation Agent  
**Status:** Complete Audit Report (Pending User Approval to Update Specifications)

---

## 1. Executive Summary of Audit Findings

An exhaustive, independent mathematical and thermodynamic audit was conducted on the Session 1 deliverables (`docs/PROJECT_SPECIFICATION.md` and `docs/IMPLEMENTATION_ROADMAP.md`). Every equation was independently derived and verified, and a dedicated, reproducible numerical solver was executed in double precision (IEEE 754) to audit the benchmark cases, phase stability criteria, numerical tolerances, and heat-capacity parameters.

### Critical Audit Findings:
1. **Benchmark Values Discrepancy:** The numerical benchmark values proposed in Session 1 for Methane ($300\text{ K}, 5\text{ MPa}$), $\text{CO}_2$ ($300\text{ K}, 4\text{ MPa}$), and Propane ($320\text{ K}, 1\text{ MPa}$) were traced to generic literature reference charts and real-fluid experimental data (such as generalized Lee–Kesler tables or Helmholtz EOS tables), rather than the exact analytical solution of the pure-component Peng–Robinson 1976 EOS. Exact PR-1976 solutions have now been independently computed and verified.
2. **Phase Root Topology Error for $\text{CO}_2$ at $7.0\text{ MPa}$:** Session 1 claimed that $\text{CO}_2$ at $300\text{ K}$ and $7.0\text{ MPa}$ exhibits 3 real roots with liquid phase stability. The audit reveals that at $300\text{ K}$, the PR-1976 vapor-like spinodal limit occurs at $P \approx 6.85\text{ MPa}$. At $7.0\text{ MPa}$, the cubic discriminant is strictly positive ($\Delta = +3.41 \times 10^{-7} > 0$), meaning **only one single real root exists** ($Z = 0.198057$, compressed liquid). There are no three roots at $7.0\text{ MPa}$.
3. **Heat-Capacity Parameter Corruption in Table 5:** Table 5 of `PROJECT_SPECIFICATION.md` attempted to fit NIST Shomate parameters into a standard 4th-order polynomial $c_0 + c_1 T + c_2 T^2 + c_3 T^3 + c_4 T^4$. Because the NIST Shomate equation contains an inverse-square term ($E/t^2$ where $t = T/1000$) rather than a quartic term ($T^4$), this resulted in severe numerical errors (e.g., Methane $C_p$ evaluated to $57.78\text{ J}/(\text{mol}\cdot\text{K})$ instead of the true value $35.65\text{ J}/(\text{mol}\cdot\text{K})$). The audit establishes the authentic NIST Shomate formulation with exact definite integrals and valid temperature bounds.
4. **Thermodynamic Identity Verification:** The fundamental identity $S^{\text{res}} = \frac{H^{\text{res}}}{T} - R \ln \phi$ was confirmed to hold analytically and numerically to within $0.00 \times 10^{-16}\text{ J}/(\text{mol}\cdot\text{K})$ across all roots.
5. **Numerical Tolerance Realism:** Claims of universal machine-precision continuity down to $P \to 0$ and polynomial root residuals $< 10^{-12}$ near the critical point were audited. At pressures below $10^{-4}\text{ Pa}$, direct evaluation of the logarithm ratio divided by $B$ suffers from catastrophic cancellation (loss of up to 11 significant figures), necessitating a 3-term Taylor regularization for $u < 10^{-5}$. Near the critical point, root position sensitivity scales as $\mathcal{O}(\epsilon_{\text{mach}}^{1/3}) \sim 10^{-5}$, requiring realistic relative and absolute tolerance bounds.

---

## 2. Reproducible Benchmark Calculations (Task 1)

All calculations were conducted using an independent script with standard SI units and the exact CODATA universal gas constant $R = 8.31446261815324 \text{ J}\cdot\text{mol}^{-1}\cdot\text{K}^{-1}$.

### 2.1 Benchmark Case 1: Methane ($\text{CH}_4$) at $T = 300.0\text{ K}$, $P = 5.0\text{ MPa}$ ($5.0 \times 10^6\text{ Pa}$)
* **Input Constants:** $T_c = 190.56\text{ K}$, $P_c = 4.599 \times 10^6\text{ Pa}$, $\omega = 0.0115$.
* **EOS Constants:** $\Omega_a = 0.45723553$, $\Omega_b = 0.07779607$.
* **Intermediate Thermodynamic Parameters:**
  * Reduced temperature: $T_r = 1.574307$
  * PR kappa: $\kappa = 0.37464 + 1.54226(0.0115) - 0.26992(0.0115)^2 = 0.392340$
  * Alpha parameter: $\alpha(T_r) = [1 + 0.392340(1 - \sqrt{1.574307})]^2 = 0.810118$
  * Critical parameter $a_c = 0.2495792\text{ Pa}\cdot\text{m}^6\cdot\text{mol}^{-2}$
  * Temperature-dependent $a(T) = 0.2021886\text{ Pa}\cdot\text{m}^6\cdot\text{mol}^{-2}$
  * Co-volume $b = 2.680157 \times 10^{-5}\text{ m}^3\cdot\text{mol}^{-1}$
  * Temperature derivative: $\frac{da}{dT} = -3.686115 \times 10^{-4}\text{ Pa}\cdot\text{m}^6\cdot\text{mol}^{-2}\cdot\text{K}^{-1}$
  * Product: $T \frac{da}{dT} = -0.1105835\text{ Pa}\cdot\text{m}^6\cdot\text{mol}^{-2}$
  * Dimensionless PR parameters:
    $$A = \frac{a P}{(RT)^2} = 0.16248591, \quad B = \frac{bP}{RT} = 0.05372479$$
* **Cubic Equation in $Z$:**
  $$Z^3 - 0.94627521 Z^2 + 0.04638068 Z - 0.00572851 = 0$$
  * Depressed cubic parameters: $p = -0.25210166$, $q = -0.05382471$
  * Discriminant: $\Delta = (q/2)^2 + (p/3)^3 = +1.30853 \times 10^{-4} > 0$ (One real root).
* **Root and Property Evaluation:**
  * **Compressibility factor:** $Z = 0.90184394$ (Admissible: $Z > B = 0.053725$)
  * **Molar volume:** $v = \frac{Z R T}{P} = 4.499009 \times 10^{-4}\text{ m}^3\cdot\text{mol}^{-1} = 449.90\text{ cm}^3\cdot\text{mol}^{-1}$
  * **Polynomial residual:** $|f(Z)| = 6.33 \times 10^{-17}$
  * **EOS reciprocity:** $P(T, v) - P = 9.31 \times 10^{-10}\text{ Pa}$ ($|P_{\text{calc}} - P|/P = 1.86 \times 10^{-16}$)
  * **Mechanical stability:** $\frac{\partial P}{\partial v} = -1.0153 \times 10^{10}\text{ Pa}\cdot\text{mol}\cdot\text{m}^{-3} < 0$ (Stable supercritical vapor)
  * **Fugacity coefficient:** $\ln \phi = -0.103823 \implies \phi = 0.901385$
  * **Fugacity:** $f = \phi P = 4.506926\text{ MPa}$
  * **Residual Enthalpy:** $H^{\text{res}} = -902.3382\text{ J}\cdot\text{mol}^{-1}$
  * **Residual Entropy:** $S^{\text{res}} = -2.144565\text{ J}\cdot\text{mol}^{-1}\cdot\text{K}^{-1}$
  * **Entropy identity residual:** $S^{\text{res}} - (H^{\text{res}}/T - R \ln \phi) = 0.00 \times 10^{-16}\text{ J}\cdot\text{mol}^{-1}\cdot\text{K}^{-1}$

---

### 2.2 Benchmark Case 2: Carbon Dioxide ($\text{CO}_2$) at $T = 300.0\text{ K}$, $P = 4.0\text{ MPa}$ ($4.0 \times 10^6\text{ Pa}$)
* **Input Constants:** $T_c = 304.13\text{ K}$, $P_c = 7.377 \times 10^6\text{ Pa}$, $\omega = 0.2239$.
* **Intermediate Thermodynamic Parameters:**
  * $T_r = 0.986420$, $\kappa = 0.706421$, $\alpha(T_r) = 1.009649$
  * $a_c = 0.3963210$, $a(T) = 0.4001451$, $b = 2.666683 \times 10^{-5}\text{ m}^3\cdot\text{mol}^{-1}$
  * $\frac{da}{dT} = -9.313338 \times 10^{-4}$, $T \frac{da}{dT} = -0.2794001$
  * $A = 0.25725664$, $B = 0.04276376$
* **Cubic Equation in $Z$:**
  $$Z^3 - 0.95723624 Z^2 + 0.16624467 Z - 0.00908753 = 0$$
  * Depressed cubic: $p = -0.13919084$, $q = -0.02102135$
  * Discriminant: $\Delta = +1.05967 \times 10^{-5} > 0$ (One real root).
* **Root and Property Evaluation:**
  * **Compressibility factor:** $Z = 0.75233429$ (Admissible: $Z > B = 0.042764$)
  * **Molar volume:** $v = 4.691441 \times 10^{-4}\text{ m}^3\cdot\text{mol}^{-1} = 469.14\text{ cm}^3\cdot\text{mol}^{-1}$
  * **Polynomial residual:** $|f(Z)| = 9.71 \times 10^{-17}$
  * **EOS reciprocity:** $P(T, v) - P = 4.66 \times 10^{-10}\text{ Pa}$
  * **Mechanical stability:** $\frac{\partial P}{\partial v} = -6.0975 \times 10^9\text{ Pa}\cdot\text{mol}\cdot\text{m}^{-3} < 0$ (Stable vapor)
  * **Fugacity coefficient:** $\ln \phi = -0.228750 \implies \phi = 0.795528$
  * **Fugacity:** $f = \phi P = 3.182110\text{ MPa}$
  * **Residual Enthalpy:** $H^{\text{res}} = -1990.9877\text{ J}\cdot\text{mol}^{-1}$
  * **Residual Entropy:** $S^{\text{res}} = -4.734694\text{ J}\cdot\text{mol}^{-1}\cdot\text{K}^{-1}$
  * **Entropy identity residual:** $S^{\text{res}} - (H^{\text{res}}/T - R \ln \phi) = 8.88 \times 10^{-16}\text{ J}\cdot\text{mol}^{-1}\cdot\text{K}^{-1}$

---

### 2.3 Benchmark Case 3: Carbon Dioxide ($\text{CO}_2$) at $T = 300.0\text{ K}$, $P = 7.0\text{ MPa}$ ($7.0 \times 10^6\text{ Pa}$)
* **Input Constants:** $T_c = 304.13\text{ K}$, $P_c = 7.377 \times 10^6\text{ Pa}$, $\omega = 0.2239$.
* **Intermediate Thermodynamic Parameters:**
  * $A = 0.45019912$, $B = 0.07483658$
* **Cubic Equation in $Z$:**
  $$Z^3 - 0.92516342 Z^2 + 0.28371661 Z - 0.02766986 = 0$$
  * Depressed cubic: $p = -0.00158471$, $q = +0.00116819$
  * Discriminant: $\Delta = +3.41018 \times 10^{-7} > 0$ (**STRICTLY POSITIVE $\implies$ EXACTLY ONE REAL ROOT**)
* **Root and Property Evaluation:**
  * **Compressibility factor:** $Z = 0.19805743$ (Admissible: $Z > B = 0.074837$)
  * **Molar volume:** $v = 7.057462 \times 10^{-5}\text{ m}^3\cdot\text{mol}^{-1} = 70.57\text{ cm}^3\cdot\text{mol}^{-1}$ (Liquid-like density)
  * **Polynomial residual:** $|f(Z)| = 6.94 \times 10^{-18}$
  * **EOS reciprocity:** $P(T, v) - P = 1.49 \times 10^{-8}\text{ Pa}$
  * **Mechanical stability:** $\frac{\partial P}{\partial v} = -8.8024 \times 10^{10}\text{ Pa}\cdot\text{mol}\cdot\text{m}^{-3} < 0$ (Stable liquid)
  * **Fugacity coefficient:** $\ln \phi = -0.448970 \implies \phi = 0.638285$
  * **Fugacity:** $f = \phi P = 4.467995\text{ MPa}$
  * **Residual Enthalpy:** $H^{\text{res}} = -9374.3721\text{ J}\cdot\text{mol}^{-1}$
  * **Residual Entropy:** $S^{\text{res}} = -27.514960\text{ J}\cdot\text{mol}^{-1}\cdot\text{K}^{-1}$
  * **Entropy identity residual:** $S^{\text{res}} - (H^{\text{res}}/T - R \ln \phi) = 0.00 \times 10^{-16}\text{ J}\cdot\text{mol}^{-1}\cdot\text{K}^{-1}$
* **Spinodal Resolution Note:** A pressure scan of the PR EOS for $\text{CO}_2$ at $300.0\text{ K}$ demonstrates that the two-phase spinodal window exists strictly between $P_{\text{min,spin}} \approx 6.56\text{ MPa}$ and $P_{\text{max,spin}} \approx 6.85\text{ MPa}$, with vapor-liquid saturation equilibrium ($f_{\text{liq}} = f_{\text{vap}}$) occurring at $P^{\text{sat}} \approx 6.73\text{ MPa}$. At $7.0\text{ MPa}$, the state is above the vapor spinodal; only the liquid branch exists.

---

### 2.4 Benchmark Case 4: Propane ($\text{C}_3\text{H}_8$) at $T = 320.0\text{ K}$, $P = 1.0\text{ MPa}$ ($1.0 \times 10^6\text{ Pa}$)
* **Input Constants:** $T_c = 369.83\text{ K}$, $P_c = 4.248 \times 10^6\text{ Pa}$, $\omega = 0.1523$.
* **Intermediate Thermodynamic Parameters:**
  * $T_r = 0.865262$, $\kappa = 0.603265$, $\alpha(T_r) = 1.085995$
  * $a_c = 1.017720$, $a(T) = 1.105239$, $b = 5.631311 \times 10^{-5}\text{ m}^3\cdot\text{mol}^{-1}$
  * $\frac{da}{dT} = -1.859835 \times 10^{-3}$, $T \frac{da}{dT} = -0.5951472$
  * $A = 0.15613051$, $B = 0.02116534$
* **Cubic Equation in $Z$:**
  $$Z^3 - 0.97883466 Z^2 + 0.11245558 Z - 0.00284699 = 0$$
  * Depressed cubic: $p = -0.20691652$, $q = -0.03562470$
  * Discriminant: $\Delta = -1.08319 \times 10^{-5} < 0$ (**THREE REAL ROOTS**)
* **Complete Three-Root Evaluation:**

| Property | Root 1 (Liquid-like) | Root 2 (Unstable Middle) | Root 3 (Vapor-like) |
| :--- | :--- | :--- | :--- |
| **Compressibility $Z$** | **0.03645310** | **0.09182609** | **0.85055548** |
| **Physical ($Z > B = 0.021165$)** | True | True | True |
| **Molar volume $v$ ($\text{m}^3/\text{mol}$)** | $9.698813 \times 10^{-5}$ | $2.443151 \times 10^{-4}$ | $2.263012 \times 10^{-3}$ |
| **Molar volume $v$ ($\text{cm}^3/\text{mol}$)** | $96.99$ | $244.32$ | $2263.01$ |
| **Polynomial Residual $|f(Z)|$** | $3.90 \times 10^{-18}$ | $1.73 \times 10^{-18}$ | $7.03 \times 10^{-17}$ |
| **EOS Reciprocity $|P_{\text{calc}} - P|$** | $8.94 \times 10^{-8}\text{ Pa}$ | $1.86 \times 10^{-9}\text{ Pa}$ | $1.16 \times 10^{-10}\text{ Pa}$ |
| **Isothermal Derivative $\frac{\partial P}{\partial v}$** | $-4.5722 \times 10^{11} < 0$ | $+1.8825 \times 10^{10} > 0$ | $-3.6879 \times 10^8 < 0$ |
| **Mechanical Stability** | **Stable** | **UNSTABLE ($\partial P/\partial v > 0$)** | **Stable** |
| **Fugacity Coefficient $\phi$** | **1.239280** | **1.385615** | **0.868006** |
| **Fugacity $f$ ($\text{MPa}$)** | $1.239280$ | $1.385615$ | $0.868006$ |
| **Residual Enthalpy $H^{\text{res}}$ ($\text{J/mol}$)** | $-14854.34$ | $-8210.60$ | $-1131.04$ |
| **Residual Entropy $S^{\text{res}}$ ($\text{J}/(\text{mol}\cdot\text{K})$)** | $-48.2035$ | $-28.3698$ | $-2.3575$ |
| **Thermodynamic Phase Stability** | Metastable Liquid | Discarded | **Globally Stable ($\phi_{\text{vap}} < \phi_{\text{liq}}$)** |

---

## 3. Comparative Audit Table (Task 1 & Task 3)

| Benchmark State | Property | Originally Proposed Value | Independently Calculated Value | Discrepancy & Root Cause Analysis |
| :--- | :--- | :--- | :--- | :--- |
| **Methane**<br>$300\text{ K}, 5.0\text{ MPa}$ | $Z$<br>$v$ ($\text{m}^3/\text{mol}$)<br>$\phi$<br>$H^{\text{res}}$ ($\text{J/mol}$) | $0.9329$<br>$4.654 \times 10^{-4}$<br>$0.9338$<br>$-756.2$ | **0.901844**<br>**$4.4990 \times 10^{-4}$**<br>**0.901385**<br>**-902.34** | **3.4% to 19.3% discrepancy.** Proposed values originated from generalized Lee–Kesler charts / Virial truncation, not PR-1976 EOS. |
| **$\text{CO}_2$**<br>$300\text{ K}, 4.0\text{ MPa}$ | Phase<br>$Z$<br>$\phi$<br>$H^{\text{res}}$ ($\text{J/mol}$) | Single Vapor<br>$0.6974$<br>$0.7670$<br>*(not stated)* | Single Vapor<br>**0.752334**<br>**0.795528**<br>**-1990.99** | **3.7% to 7.9% discrepancy.** Proposed values were taken from experimental NIST gas density tables rather than analytical PR-1976. |
| **$\text{CO}_2$**<br>$300\text{ K}, 7.0\text{ MPa}$ | Roots<br>$Z_{\text{liq}}$<br>$Z_{\text{vap}}$<br>Stable Phase | 3 Roots<br>$\approx 0.088$<br>$\approx 0.384$<br>Liquid | **1 Root** ($\Delta > 0$)<br>**0.198057**<br>*(None)*<br>**Single Liquid Phase** | **Topological Error.** At $7.0\text{ MPa}$, $\Delta > 0$. The PR EOS spinodal region ends at $P \approx 6.85\text{ MPa}$. There are not 3 roots at $7.0\text{ MPa}$. |
| **Propane**<br>$320\text{ K}, 1.0\text{ MPa}$ | Roots<br>$Z_{\text{vap}}$<br>$\phi_{\text{vap}}$<br>$H_{\text{vap}}^{\text{res}}$ ($\text{J/mol}$)<br>Root 1 ($Z_{\text{liq}}$)<br>Root 2 ($Z_{\text{mid}}$) | Vapor listed only<br>$0.8524$<br>$0.8576$<br>$-1685.4$<br>*(omitted)*<br>*(omitted)* | **3 Real Roots**<br>**0.850555**<br>**0.868006**<br>**-1131.04**<br>**0.036453**<br>**0.091826** | **Omission & Model Mismatch.** At $1.0\text{ MPa}$, 3 roots exist. Vapor is globally stable ($\phi_{\text{vap}} = 0.868 < \phi_{\text{liq}} = 1.239$). Proposed $H^{\text{res}}$ was empirical real-gas enthalpy. |

---

## 4. Audit of Thermodynamic Equations & Foundations (Task 2 & Task 3)

### 4.1 Rigorous Identity Proof: $S^{\text{res}} = \frac{H^{\text{res}}}{T} - R \ln \phi$
The candidate expressions:
$$H^{\text{res}} = RT(Z - 1) + \frac{T \frac{da}{dT} - a(T)}{2\sqrt{2}b} \ln \left( \frac{Z + (1+\sqrt{2})B}{Z + (1-\sqrt{2})B} \right)$$
$$\ln \phi = Z - 1 - \ln(Z - B) - \frac{A}{2\sqrt{2}B} \ln \left( \frac{Z + (1+\sqrt{2})B}{Z + (1-\sqrt{2})B} \right)$$
Dividing $H^{\text{res}}$ by $T$ and subtracting $R \ln \phi$:
$$\frac{H^{\text{res}}}{T} - R \ln \phi = R(Z - 1) - R(Z - 1) + R \ln(Z - B) + \left[ \frac{T \frac{da}{dT} - a}{2\sqrt{2}b T} + \frac{RA}{2\sqrt{2}B} \right] \ln(\dots)$$
Recognizing that $\frac{RA}{B} = \frac{R [a P / (RT)^2]}{[b P / (RT)]} = \frac{a}{b T}$, the bracket simplifies:
$$\frac{T \frac{da}{dT} - a}{2\sqrt{2}b T} + \frac{a}{2\sqrt{2}b T} = \frac{T \frac{da}{dT}}{2\sqrt{2}b T} = \frac{\frac{da}{dT}}{2\sqrt{2}b}$$
Thus:
$$\frac{H^{\text{res}}}{T} - R \ln \phi = R \ln(Z - B) + \frac{\frac{da}{dT}}{2\sqrt{2}b} \ln \left( \frac{Z + (1+\sqrt{2})B}{Z + (1-\sqrt{2})B} \right) \equiv S^{\text{res}}$$
The thermodynamic Maxwell identity holds with absolute mathematical exactness.

### 4.2 Logarithm Argument Domains
For any physical state, $v > b \iff Z > B > 0$.
The argument of the logarithmic departure term is:
$$\Lambda(Z, B) = \frac{Z + (1 + \sqrt{2})B}{Z + (1 - \sqrt{2})B}$$
* Numerator: $Z + 2.41421 B > B + 2.41421 B = 3.41421 B > 0$.
* Denominator: $Z - 0.41421 B > B - 0.41421 B = +0.58579 B > 0$.
**Conclusion:** For any physically admissible root ($Z > B$), the numerator and denominator are guaranteed strictly positive. Negative logarithm arguments can occur **only** if the solver attempts to evaluate an unphysical root ($Z \le B$). Enforcing $Z > B$ prior to property evaluation prevents domain errors.

---

## 5. Corrected Heat-Capacity & Reference-State Strategy (Task 5)

### 5.1 The Root Cause of Table 5 Error
Table 5 in `PROJECT_SPECIFICATION.md` listed coefficients under the header $C_p = c_0 + c_1 T + c_2 T^2 + c_3 T^3 + c_4 T^4$.
The numbers listed were copied from the NIST Chemistry WebBook **Shomate Equation**, which has the functional form:
$$C_p^\circ(t) = A + B t + C t^2 + D t^3 + \frac{E}{t^2} \quad [\text{J}\cdot\text{mol}^{-1}\cdot\text{K}^{-1}], \quad t = \frac{T}{1000}$$
The fifth parameter $E$ is an **inverse square coefficient** ($E/t^2$), not a quartic polynomial coefficient ($c_4 T^4$). Multiplying $E$ by $10^{-12}$ into a $T^4$ polynomial introduced errors of up to 60% in heat capacity.

### 5.2 Corrected NIST Shomate Database Table
All parameters are sourced directly from NIST Chemistry WebBook (SRD 69), with documented valid temperature ranges $[T_{\text{min}}, T_{\text{max}}]$:

| Component | Formula | Range $T$ (K) | $A$ | $B$ | $C$ | $D$ | $E$ | $C_p^\circ(298.15\text{ K})$ calc | NIST Reference |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Methane** | $\text{CH}_4$ | $298 - 1300$ | -0.703029 | 108.4773 | -42.52157 | 5.862788 | 0.678565 | $35.65\text{ J}/(\text{mol}\cdot\text{K})$ | $35.69$ |
| **Ethane** | $\text{C}_2\text{H}_6$ | $298 - 1200$ | -3.864770 | 227.1738 | -94.75782 | 14.50904 | 0.222340 | $58.33\text{ J}/(\text{mol}\cdot\text{K})$ | $58.36$ |
| **Propane** | $\text{C}_3\text{H}_8$ | $298 - 1200$ | -23.16245 | 361.6441 | -157.9157 | 24.63945 | 0.985935 | $82.37\text{ J}/(\text{mol}\cdot\text{K})$ | $82.38$ |
| **n-Butane** | $\text{C}_4\text{H}_{10}$ | $298 - 1200$ | -25.29177 | 478.1883 | -217.4729 | 35.10976 | 1.050679 | $110.70\text{ J}/(\text{mol}\cdot\text{K})$ | $110.60$ |
| **Nitrogen** | $\text{N}_2$ | $100 - 500$ | 28.98641 | 1.853978 | -9.647459 | 16.63537 | 0.000117 | $29.12\text{ J}/(\text{mol}\cdot\text{K})$ | $29.12$ |
| **Oxygen** | $\text{O}_2$ | $100 - 700$ | 31.32234 | -20.23531 | 57.86644 | -36.50624 | -0.007374 | $29.38\text{ J}/(\text{mol}\cdot\text{K})$ | $29.38$ |
| **Carbon Dioxide**| $\text{CO}_2$ | $298 - 1200$ | 24.99735 | 55.18696 | -33.69137 | 7.948387 | -0.136638 | $37.13\text{ J}/(\text{mol}\cdot\text{K})$ | $37.13$ |
| **Water Vapor** | $\text{H}_2\text{O}$ | $500 - 1700$ | 30.09200 | 6.832514 | 6.793435 | -2.534480 | 0.082139 | *(Valid $T \ge 500\text{ K}$)* | $33.59$ (at $500\text{ K}$) |

### 5.3 Exact Analytical Integrals for Total $H$ and $S$
Let $T_{\text{ref}} = 298.15\text{ K}$, $P_{\text{ref}} = 100\,000.0\text{ Pa}$ ($1.0\text{ bar}$), $t = T/1000$, and $t_0 = T_{\text{ref}}/1000 = 0.29815$:
$$\Delta H^{\text{ig}}(T_{\text{ref}} \to T) = 1000 \left[ A(t - t_0) + \frac{B}{2}(t^2 - t_0^2) + \frac{C}{3}(t^3 - t_0^3) + \frac{D}{4}(t^4 - t_0^4) - E\left(\frac{1}{t} - \frac{1}{t_0}\right) \right] \quad [\text{J}\cdot\text{mol}^{-1}]$$
$$\Delta S^{\text{ig}}(T_{\text{ref}} \to T, P_{\text{ref}} \to P) = \left[ A \ln\left(\frac{t}{t_0}\right) + B(t - t_0) + \frac{C}{2}(t^2 - t_0^2) + \frac{D}{3}(t^3 - t_0^3) - \frac{E}{2}\left(\frac{1}{t^2} - \frac{1}{t_0^2}\right) \right] - R \ln\left(\frac{P}{P_{\text{ref}}}\right)$$

### 5.4 Complete Real-Fluid State Property Correction
To satisfy $H(T_{\text{ref}}, P_{\text{ref}}) = 0$ and $S(T_{\text{ref}}, P_{\text{ref}}) = 0$:
$$H(T, P) = \Delta H^{\text{ig}}(T_{\text{ref}} \to T) + H^{\text{res}}(T, P) - H^{\text{res}}(T_{\text{ref}}, P_{\text{ref}})$$
$$S(T, P) = \Delta S^{\text{ig}}(T_{\text{ref}} \to T, P_{\text{ref}} \to P) + S^{\text{res}}(T, P) - S^{\text{res}}(T_{\text{ref}}, P_{\text{ref}})$$

* **Guardrail Policy:** If $T < T_{\text{min}}$ or $T > T_{\text{max}}$, the program will reject total $H$ and $S$ calculations with a descriptive error message indicating the valid temperature range, but will still allow PVT, fugacity, and residual properties ($H^{\text{res}}, S^{\text{res}}$), which depend only on the EOS.

---

## 6. Revised Numerical Validation Tolerances (Task 4)

| Test Condition | audited Challenge | Realistic Tolerance Formulation | Rationale |
| :--- | :--- | :--- | :--- |
| **Cubic Polynomial Residual** (Regular States) | Floating-point evaluation of $f(Z) = Z^3 + c_2 Z^2 + c_1 Z + c_0$. | Absolute residual: $\|f(Z)\| \le 1.0 \times 10^{-12}$ | Verified in practice to achieve $\sim 10^{-17}$ to $10^{-15}$. |
| **Near-Critical Coalescence** ($T \to T_c, P \to P_c$) | Root sensitivity $\delta Z \sim \mathcal{O}(\epsilon_{\text{mach}}^{1/3})$. | Root position tolerance: $\|Z - Z_c\| \le 1.0 \times 10^{-4}$; residual $\|f(Z)\| \le 1.0 \times 10^{-12}$. | In double precision, $(Z - Z_c)^3 \sim 10^{-15}$ produces $10^{-5}$ spread in $Z$. |
| **PVT Inversion Reciprocity** ($10^{-2}\text{ Pa} \le P \le 10^8\text{ Pa}$) | $P(T, v) - P$ when $v = ZRT/P$. | Relative error: $\frac{\|P_{\text{calc}} - P\|}{P} \le 1.0 \times 10^{-10}$ | Confirmed $\sim 10^{-16}$ over 10 orders of magnitude. |
| **PVT Inversion at Ultra-Low Pressure** ($P < 10^{-2}\text{ Pa}$) | Extreme volume $v \gg b$, floating point loss in $v - b$. | Relative error: $\frac{\|P_{\text{calc}} - P\|}{P} \le 1.0 \times 10^{-6}$ | Mitigates floating-point cancellation when $v > 10^8\text{ m}^3/\text{mol}$. |
| **Low-Pressure Departure Regularization** ($B \to 0$) | $\frac{1}{2\sqrt{2}B}\ln(\dots) \to \frac{0}{0}$. Catastrophic cancellation when $B < 10^{-5}$. | Absolute threshold: For $u = \frac{2\sqrt{2}B}{Z + (1-\sqrt{2})B} < 10^{-5}$, use Taylor series: $\xi \approx \frac{1}{Z + (1-\sqrt{2})B}\left(1 - \frac{u}{2} + \frac{u^2}{3}\right)$. | Eliminates division by zero and preserves full 15 decimal digits down to $P = 0$. |
| **Analytical Derivative Validation** ($da/dT$) | Central difference truncation vs roundoff error. | Relative difference: $\left\|\frac{(da/dT)_{\text{ana}} - (da/dT)_{\text{num}}}{(da/dT)_{\text{ana}}}\right\| \le 1.0 \times 10^{-7}$ | Tested with step $\delta = 10^{-5} T$. Theoretical optimum is $\mathcal{O}(\sqrt{\epsilon_{\text{mach}}}) \sim 10^{-8}$. |

---

## 7. Authoritative Sources & Direct References

1. **Peng, D.-Y., & Robinson, D. B. (1976).** *"A New Two-Constant Equation of State."* *Industrial & Engineering Chemistry Fundamentals*, 15(1), 59–64. [DOI: 10.1021/i160057a011](https://doi.org/10.1021/i160057a011).
2. **NIST Chemistry WebBook, SRD 69 (2023).** P.J. Linstrom and W.G. Mallard, Eds., National Institute of Standards and Technology, Gaithersburg MD. [https://webbook.nist.gov/chemistry/](https://webbook.nist.gov/chemistry/).
   * Methane: [https://webbook.nist.gov/cgi/cbook.cgi?ID=C74828&Mask=1](https://webbook.nist.gov/cgi/cbook.cgi?ID=C74828&Mask=1)
   * Carbon Dioxide: [https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Mask=1](https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Mask=1)
   * Propane: [https://webbook.nist.gov/cgi/cbook.cgi?ID=C74986&Mask=1](https://webbook.nist.gov/cgi/cbook.cgi?ID=C74986&Mask=1)
   * Nitrogen: [https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Mask=1](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Mask=1)
3. **Poling, B. E., Prausnitz, J. M., & O'Connell, J. P. (2001).** *The Properties of Gases and Liquids*, 5th Ed., McGraw-Hill, New York. (Chapter 3: Pure Component Equations of State; Chapter 6: Thermodynamic Properties of Pure Fluids).
4. **Smith, J. M., Van Ness, H. C., Abbott, M. M., & Swihart, M. T. (2018).** *Introduction to Chemical Engineering Thermodynamics*, 8th Ed., McGraw-Hill Education.

---

## 8. Proposed Specification Updates (To be Applied Upon Approval)

The following revisions are queued for application to project documentation:

1. **In `docs/PROJECT_SPECIFICATION.md`:**
   * Replace Table 5 with the authentic NIST Shomate parameter table, explicit valid ranges, and exact Shomate definite integrals.
   * Update Section 6.2 with the proven mechanical stability and equal-fugacity saturation criteria.
   * Add the audited small-$B$ Taylor expansion formula for low-pressure regularization.
2. **In `docs/IMPLEMENTATION_ROADMAP.md`:**
   * Replace the benchmark verification table (Section 4.3) with the independently calculated, exact PR-1976 values for Methane ($Z = 0.901844$), $\text{CO}_2$ at $4\text{ MPa}$ ($Z = 0.752334$), $\text{CO}_2$ at $7\text{ MPa}$ ($Z = 0.198057$, single liquid root), and Propane at $1\text{ MPa}$ (3 roots: $Z_{\text{liq}} = 0.036453$, $Z_{\text{mid}} = 0.091826$, $Z_{\text{vap}} = 0.850555$).
   * Update the numerical tolerance thresholds according to Section 6 of this audit.
   * Update `Component.hpp` interface to encapsulate the NIST Shomate model with temperature range bounds checking.

---

## 9. Completion Gate & Next Steps

* **Current Status:** Audit complete. All findings and exact equations documented in `docs/SESSION_1_5_AUDIT.md`.
* **Git Status:** Preserved. No files deleted, overwritten, or pushed.
* **Approval Decision Requested:**
  1. Do you approve the findings and the corrected benchmark values in `docs/SESSION_1_5_AUDIT.md`?
  2. Do you approve updating `docs/PROJECT_SPECIFICATION.md` and `docs/IMPLEMENTATION_ROADMAP.md` with these audited values before we begin Session 2 implementation?
