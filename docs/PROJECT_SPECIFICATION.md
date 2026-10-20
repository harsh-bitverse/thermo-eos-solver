# CHL2002 Thermodynamics Project: System Specification & Scientific Foundations

**Course:** CHL2002 Thermodynamics  
**Project Title:** Pure-Component Real Gas PVT, Departure Properties, and Phase Fugacity via Peng–Robinson Equation of State  
**Language & Standard:** C++17 / C++20 (Standard Library only)  
**Build System:** CMake (>= 3.20)  
**Submission Due Date:** November 10, 2026  
**Document Status:** Approved Specification (Session 1)

---

## 1. Executive Summary and Scope

This software application provides a scientifically rigorous, numerically robust C++ implementation of the pure-component Peng–Robinson (PR) cubic equation of state (1976). Designed for chemical engineering thermodynamics, the system evaluates:
1. **PVT Behavior:** Direct calculation of molar volume $v$ and compressibility factor $Z$ given $(T, P)$, and pressure $P$ given $(T, v)$.
2. **Phase Root Identification:** Complete analytical resolution of cubic roots, physical admissibility filtering ($Z > B \iff v > b$), and phase stability determination via chemical potential / fugacity minimization.
3. **Residual (Departure) Properties:** Exact analytical computation of residual enthalpy $H^{\text{res}}(T, P)$ and residual entropy $S^{\text{res}}(T, P)$.
4. **Total State Properties:** Rigorous integration of ideal-gas heat capacity $C_p^{\text{ig}}(T)$ relative to a standard thermodynamic reference state ($T_{\text{ref}} = 298.15\text{ K}, P_{\text{ref}} = 10^5\text{ Pa}$).
5. **Pure-Component Fugacity:** Calculation of fugacity coefficient $\phi(T, P)$ and fugacity $f(T, P)$ with analytical asymptotic regularization in the ideal-gas limit ($P \to 0$).
6. **Curated & Extensible Database:** Traceable critical properties and heat-capacity polynomial coefficients for common industrial gases, plus support for custom user-defined components.

---

## 2. Requirements Mapping & Traceability Matrix

| Requirement ID | Specification Item | Mathematical / Functional Requirement | Acceptance Criteria |
| :--- | :--- | :--- | :--- |
| **REQ-PVT-01** | Forward PVT ($(T, P) \to Z, v$) | Solve cubic in $Z$, compute $v = Z R T / P$. | All physical roots identified; $v > b$. |
| **REQ-PVT-02** | Inverse PVT ($(T, v) \to P$) | Direct evaluation of $P(T, v) = \frac{RT}{v-b} - \frac{a(T)}{v^2 + 2bv - b^2}$. | Matches forward calculation within $10^{-12}\text{ Pa}$. |
| **REQ-ROOT-01** | Cubic Root Solution | Solve $Z^3 - (1-B)Z^2 + (A - 3B^2 - 2B)Z - (AB - B^2 - B^3) = 0$. | Analytic Cardano/trigonometric method; residual $< 10^{-12}$. |
| **REQ-ROOT-02** | Root Classification & Stability | Identify vapor ($Z_{\max}$), liquid ($Z_{\min}$), and mechanically unstable intermediate root. Compare $G^{\text{res}}$ (or $\ln \phi$) to establish thermodynamically stable phase. | Unphysical roots ($Z \le B$) discarded; stable phase flagged. |
| **REQ-DEP-01** | Residual Enthalpy $H^{\text{res}}$ | Exact analytical derivation with $da/dT$. | $\lim_{P\to 0} H^{\text{res}} = 0$; finite-difference check on $da/dT < 10^{-7}$. |
| **REQ-DEP-02** | Residual Entropy $S^{\text{res}}$ | Exact analytical derivation with $da/dT$ and $R\ln(Z-B)$. | $\lim_{P\to 0} S^{\text{res}} = 0$. |
| **REQ-FUG-01** | Fugacity & Fugacity Coefficient | $\ln \phi = Z - 1 - \ln(Z-B) - \frac{A}{2\sqrt{2}B}\ln[\dots]$; $f = \phi P$. | $\lim_{P\to 0} \phi = 1.0$; $f \to P$. Numerically regularized for $B \to 0$. |
| **REQ-CP-01** | Ideal Gas Integration | Fourth-order polynomial $C_p^{\text{ig}}(T) = \sum_{i=0}^4 c_i T^i$. Analytical integrals for $\Delta H^{\text{ig}}$ and $\Delta S^{\text{ig}}$. | Explicit reference state ($298.15\text{ K}, 1\text{ bar}$) documented. |
| **REQ-DATA-01** | Component Catalog | Built-in data for $\ge 6$ standard gases with traceable citations. | Values sourced from NIST WebBook and Poling et al. |
| **REQ-DATA-02** | Custom Component Entry | Runtime user entry of $T_c, P_c, \omega$, and optional $C_p^{\text{ig}}$ coefficients. | Input validation ($T_c > 0, P_c > 0, \omega > -1$). |
| **REQ-ARCH-01** | Modular C++ Architecture | Separation of thermo core, numerics, database, and CLI. | Pure C++17, zero external dependencies, standard CMake. |
| **REQ-TEST-01** | Automated Test Suite | Unit tests for cubic solver, limits, reciprocity, and literature benchmarks. | All tests pass reproducibly on Windows MinGW/GCC and MSVC. |

---

## 3. Verified Thermodynamic Formulations

### 3.1 Peng–Robinson Equation of State (1976)
The pure-component Peng–Robinson EOS expresses pressure $P$ as:
$$P = \frac{RT}{v - b} - \frac{a(T)}{v(v + b) + b(v - b)} = \frac{RT}{v - b} - \frac{a(T)}{v^2 + 2bv - b^2}$$

where:
- $R = 8.31446261815324 \text{ J}\cdot\text{mol}^{-1}\cdot\text{K}^{-1}$ (CODATA 2018/2022 exact value).
- $T$ is absolute temperature ($\text{K}$).
- $v$ is molar volume ($\text{m}^3\cdot\text{mol}^{-1}$).
- $P$ is pressure ($\text{Pa}$).

The parameters $a(T)$ and $b$ are determined from critical constants ($T_c, P_c$) and acentric factor $\omega$:
$$a_c = \Omega_a \frac{R^2 T_c^2}{P_c}, \quad \Omega_a = \frac{0.4572355289...}{\approx 0.45724}$$
$$b = \Omega_b \frac{R T_c}{P_c}, \quad \Omega_b = \frac{0.0777960739...}{\approx 0.07780}$$
$$a(T) = a_c \cdot \alpha(T_r, \omega)$$
$$T_r = \frac{T}{T_c}$$
$$\alpha(T_r, \omega) = \left[ 1 + \kappa \left( 1 - \sqrt{T_r} \right) \right]^2$$
$$\kappa = 0.37464 + 1.54226\,\omega - 0.26992\,\omega^2$$

*(Note on PR-1978: For heavier hydrocarbons where $\omega > 0.49$, Peng and Robinson (1978) introduced $\kappa = 0.379642 + 1.48503\omega - 0.164423\omega^2 + 0.016666\omega^3$. In this project, the standard 1976 definition is implemented as specified by the course syllabus).*

### 3.2 Temperature Derivative of $a(T)$
The temperature derivative $\frac{da}{dT}$ is required for departure enthalpy and entropy:
$$\alpha(T) = \left[ 1 + \kappa \left( 1 - \sqrt{\frac{T}{T_c}} \right) \right]^2$$
$$\frac{d\alpha}{dT} = 2 \left[ 1 + \kappa \left( 1 - \sqrt{T_r} \right) \right] \left( -\frac{\kappa}{2 \sqrt{T \cdot T_c}} \right) = -\frac{\kappa \sqrt{\alpha(T)}}{\sqrt{T \cdot T_c}}$$
$$\frac{da}{dT} = a_c \frac{d\alpha}{dT} = - \frac{a_c \kappa \sqrt{\alpha(T)}}{\sqrt{T \cdot T_c}} = - \frac{a_c \kappa \left[ 1 + \kappa \left( 1 - \sqrt{T_r} \right) \right]}{\sqrt{T \cdot T_c}}$$

In product form:
$$T \frac{da}{dT} = - a_c \kappa \sqrt{T_r} \left[ 1 + \kappa \left( 1 - \sqrt{T_r} \right) \right] = - a_c \kappa \sqrt{T_r \alpha(T)}$$

### 3.3 Compressibility Factor Polynomial
Substituting $v = \frac{Z R T}{P}$ into the PR EOS yields the dimensionless cubic in $Z$:
$$Z^3 - (1 - B)Z^2 + (A - 3B^2 - 2B)Z - (AB - B^2 - B^3) = 0$$

where dimensionless parameters $A$ and $B$ are:
$$A = \frac{a(T) P}{(R T)^2}, \quad B = \frac{b P}{R T}$$

Polynomial coefficients:
$$c_3 = 1$$
$$c_2 = -(1 - B) = B - 1$$
$$c_1 = A - 3B^2 - 2B$$
$$c_0 = -(AB - B^2 - B^3)$$

### 3.4 Rigorous Derivation of Residual (Departure) Properties
A residual property $M^{\text{res}}(T, P)$ is defined as the difference between the actual property of a real fluid and that of an ideal gas at the same temperature and pressure:
$$M^{\text{res}}(T, P) \equiv M(T, P) - M^{\text{ig}}(T, P)$$

#### 3.4.1 Residual Gibbs Energy and Fugacity Coefficient
From fundamental thermodynamics:
$$\frac{G^{\text{res}}(T, P)}{RT} = \ln \phi = \int_0^P (Z - 1) \frac{dP'}{P'} = Z - 1 - \ln Z - \frac{1}{RT} \int_\infty^v \left( P - \frac{RT}{v'} \right) dv'$$

Evaluating the volume integral for the PR EOS:
$$\int_\infty^v \left( \frac{RT}{v' - b} - \frac{RT}{v'} \right) dv' = RT \ln \left( \frac{v - b}{v} \right) = RT \ln \left( 1 - \frac{b}{v} \right)$$
$$\int_\infty^v \frac{a(T)}{v'^2 + 2bv' - b^2} dv' = \frac{a(T)}{2\sqrt{2}b} \left[ \ln \left( \frac{v + (1-\sqrt{2})b}{v + (1+\sqrt{2})b} \right) \right]_\infty^v = \frac{a(T)}{2\sqrt{2}b} \ln \left( \frac{v + (1-\sqrt{2})b}{v + (1+\sqrt{2})b} \right)$$

Combining terms and noting $\ln Z + \ln\left(\frac{v-b}{v}\right) = \ln\left(\frac{P(v-b)}{RT}\right) = \ln(Z - B)$:
$$\ln \phi = Z - 1 - \ln(Z - B) - \frac{A}{2\sqrt{2} B} \ln \left( \frac{Z + (1 + \sqrt{2})B}{Z + (1 - \sqrt{2})B} \right)$$

Fugacity $f$ is:
$$f = \phi \cdot P \quad [\text{Pa}]$$

#### 3.4.2 Residual Enthalpy
From the Gibbs–Helmholtz relation $\frac{H^{\text{res}}}{RT^2} = - \left( \frac{\partial (G^{\text{res}}/RT)}{\partial T} \right)_P$, or via volume integration:
$$H^{\text{res}}(T, P) = \int_\infty^v \left[ T \left( \frac{\partial P}{\partial T} \right)_v - P \right] dv' + RT(Z - 1)$$

For PR EOS:
$$T \left( \frac{\partial P}{\partial T} \right)_v - P = \frac{a(T) - T \frac{da}{dT}}{v^2 + 2bv - b^2}$$
$$\int_\infty^v \frac{a - T \frac{da}{dT}}{v'^2 + 2bv' - b^2} dv' = \frac{T \frac{da}{dT} - a}{2\sqrt{2}b} \ln \left( \frac{Z + (1 + \sqrt{2})B}{Z + (1 - \sqrt{2})B} \right)$$

Therefore:
$$H^{\text{res}}(T, P) = RT(Z - 1) + \frac{T \frac{da}{dT} - a(T)}{2\sqrt{2}b} \ln \left( \frac{Z + (1 + \sqrt{2})B}{Z + (1 - \sqrt{2})B} \right) \quad [\text{J}\cdot\text{mol}^{-1}]$$

#### 3.4.3 Residual Entropy
Using $S^{\text{res}} = \frac{H^{\text{res}} - G^{\text{res}}}{T} = \frac{H^{\text{res}}}{T} - R \ln \phi$:
$$\frac{H^{\text{res}}}{T} = R(Z - 1) + \frac{T \frac{da}{dT} - a}{2\sqrt{2}b T} \ln \left( \frac{Z + (1 + \sqrt{2})B}{Z + (1 - \sqrt{2})B} \right)$$
$$-R \ln \phi = -R(Z - 1) + R \ln(Z - B) + \frac{R A}{2\sqrt{2} B} \ln \left( \frac{Z + (1 + \sqrt{2})B}{Z + (1 - \sqrt{2})B} \right)$$

Since $\frac{R A}{B} = \frac{a(T)}{b T}$, the bracketed coefficient simplifies identically:
$$\frac{T \frac{da}{dT} - a}{2\sqrt{2}b T} + \frac{a}{2\sqrt{2}b T} = \frac{T \frac{da}{dT}}{2\sqrt{2}b T} = \frac{\frac{da}{dT}}{2\sqrt{2}b}$$

Therefore:
$$S^{\text{res}}(T, P) = R \ln(Z - B) + \frac{\frac{da}{dT}}{2\sqrt{2}b} \ln \left( \frac{Z + (1 + \sqrt{2})B}{Z + (1 - \sqrt{2})B} \right) \quad [\text{J}\cdot\text{mol}^{-1}\cdot\text{K}^{-1}]$$

---

## 4. Ideal-Gas Heat Capacity & Reference-State Convention

### 4.1 Heat Capacity Model
Ideal-gas heat capacity $C_p^{\text{ig}}(T)$ is parameterized as a 4th-order polynomial:
$$C_p^{\text{ig}}(T) = c_0 + c_1 T + c_2 T^2 + c_3 T^3 + c_4 T^4 \quad [\text{J}\cdot\text{mol}^{-1}\cdot\text{K}^{-1}]$$

### 4.2 Exact Integrals for Enthalpy and Entropy
Let $T_{\text{ref}} = 298.15\text{ K}$ and $P_{\text{ref}} = 100\,000\text{ Pa}$ (1.0 bar).
$$\Delta H^{\text{ig}}(T_{\text{ref}} \to T) = \int_{T_{\text{ref}}}^T C_p^{\text{ig}}(T') dT' = \sum_{k=0}^4 \frac{c_k}{k+1} \left( T^{k+1} - T_{\text{ref}}^{k+1} \right)$$
$$\Delta S^{\text{ig}}(T_{\text{ref}} \to T, P_{\text{ref}} \to P) = \int_{T_{\text{ref}}}^T \frac{C_p^{\text{ig}}(T')}{T'} dT' - R \ln \left( \frac{P}{P_{\text{ref}}} \right)$$
$$\int_{T_{\text{ref}}}^T \frac{C_p^{\text{ig}}(T')}{T'} dT' = c_0 \ln \left( \frac{T}{T_{\text{ref}}} \right) + \sum_{k=1}^4 \frac{c_k}{k} \left( T^k - T_{\text{ref}}^k \right)$$

### 4.3 Total Property Reporting Policy
To prevent ambiguity in technical reports and exams:
1. **Residual Properties:** Always report $H^{\text{res}}(T, P)$ and $S^{\text{res}}(T, P)$ directly. These are absolute departures that depend only on the EOS and state $(T, P)$, with zero reference-state ambiguity.
2. **Total Relative Properties:** Report $H(T, P) - H^{\circ}(T_{\text{ref}}, P_{\text{ref}})$ and $S(T, P) - S^{\circ}(T_{\text{ref}}, P_{\text{ref}})$, defined explicitly as:
   $$H(T, P) = \Delta H^{\text{ig}}(T_{\text{ref}} \to T) + H^{\text{res}}(T, P) - H^{\text{res}}(T_{\text{ref}}, P_{\text{ref}})$$
   $$S(T, P) = \Delta S^{\text{ig}}(T_{\text{ref}} \to T, P_{\text{ref}} \to P) + S^{\text{res}}(T, P) - S^{\text{res}}(T_{\text{ref}}, P_{\text{ref}})$$
   under the reference convention $H(T_{\text{ref}}, P_{\text{ref}}) = 0$ and $S(T_{\text{ref}}, P_{\text{ref}}) = 0$.

---

## 5. Curated Component Database

Parameters compiled from *Poling, Prausnitz, & O'Connell (2001), The Properties of Gases and Liquids (5th Ed.)* and *NIST Chemistry WebBook (SRD 69)*:

| Component | Formula | $M$ (g/mol) | $T_c$ (K) | $P_c$ (bar) | $P_c$ (MPa) | $\omega$ | $c_0$ | $c_1 \times 10^2$ | $c_2 \times 10^5$ | $c_3 \times 10^9$ | $c_4 \times 10^{12}$ |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Methane** | $\text{CH}_4$ | 16.043 | 190.56 | 45.99 | 4.599 | 0.0115 | 34.398 | -5.586 | 56.40 | -410.0 | 97.6 |
| **Ethane** | $\text{C}_2\text{H}_6$ | 30.070 | 305.32 | 48.72 | 4.872 | 0.0995 | 28.140 | 8.644 | 24.71 | -220.0 | 56.8 |
| **Propane** | $\text{C}_3\text{H}_8$ | 44.096 | 369.83 | 42.48 | 4.248 | 0.1523 | 31.920 | 14.36 | 31.80 | -330.0 | 88.5 |
| **n-Butane** | $\text{C}_4\text{H}_{10}$ | 58.122 | 425.12 | 37.96 | 3.796 | 0.2002 | 36.160 | 21.65 | 39.10 | -420.0 | 114.0 |
| **Nitrogen** | $\text{N}_2$ | 28.013 | 126.20 | 33.98 | 3.398 | 0.0372 | 29.105 | -0.192 | 0.400 | -0.871 | 0.0 |
| **Oxygen** | $\text{O}_2$ | 31.999 | 154.58 | 50.43 | 5.043 | 0.0222 | 29.530 | -0.889 | 3.808 | -3.262 | 0.886 |
| **Carbon Dioxide**| $\text{CO}_2$ | 44.010 | 304.13 | 73.77 | 7.377 | 0.2239 | 24.997 | 5.519 | -3.369 | 7.948 | -0.137 |
| **Water** | $\text{H}_2\text{O}$ | 18.015 | 647.10 | 220.64| 22.064 | 0.3443 | 32.240 | 0.192 | 1.055 | -0.360 | 0.0 |

*Internal conversion:* Pressures are stored and computed in SI units ($\text{Pa} = 10^5 \times \text{bar} = 10^6 \times \text{MPa}$).

---

## 6. Numerical Methods & Algorithm Design

### 6.1 Cubic Equation Solver
Given $Z^3 + c_2 Z^2 + c_1 Z + c_0 = 0$, substitute $Z = x - \frac{c_2}{3}$ to obtain the depressed cubic:
$$x^3 + px + q = 0$$
where:
$$p = c_1 - \frac{c_2^2}{3}, \quad q = c_0 - \frac{c_1 c_2}{3} + \frac{2 c_2^3}{27}$$

The discriminant is:
$$\Delta = \left( \frac{q}{2} \right)^2 + \left( \frac{p}{3} \right)^3$$

#### Case 1: $\Delta < 0$ (Three distinct real roots)
Occurs in the two-phase region. To prevent round-off error from complex cube roots, use Vieta's trigonometric formulation:
$$\theta = \frac{1}{3} \arccos \left( \frac{-q/2}{\sqrt{-(p/3)^3}} \right)$$
$$x_k = 2 \sqrt{-\frac{p}{3}} \cos \left( \theta + \frac{2k\pi}{3} \right), \quad k \in \{0, 1, 2\}$$
$$Z_k = x_k - \frac{c_2}{3}$$

#### Case 2: $\Delta \ge 0$ (One real root, or multiple equal roots)
Cardano's analytical formulation:
$$u = \sqrt[3]{ -\frac{q}{2} + \sqrt{\Delta} }, \quad v = \sqrt[3]{ -\frac{q}{2} - \sqrt{\Delta} }$$
$$x_1 = u + v \implies Z_1 = x_1 - \frac{c_2}{3}$$

### 6.2 Root Filtering & Physical Root Selection
1. **Physical Admissibility:** Molar volume must satisfy $v > b$, which requires $Z > B$. Any root with $Z \le B$ is non-physical and immediately rejected.
2. **Single Root ($T > T_c$ or outside two-phase envelope):** A single physical root exists, corresponding to supercritical fluid or single-phase gas/liquid.
3. **Three Real Roots ($T < T_c$ within spinodal limits):**
   - Largest root: $Z_{\text{vap}} = \max(Z_k)$ (Vapor-like phase).
   - Smallest root: $Z_{\text{liq}} = \min(Z_k)$ (Liquid-like phase).
   - Intermediate root: Mechanically unstable ($\partial P / \partial v > 0$), physically unobservable.
4. **Thermodynamic Phase Stability:**
   At $(T, P)$, the thermodynamically stable phase is the one that minimizes molar Gibbs energy $G(T, P)$. Since $G^{\text{ig}}$ is identical for all roots at $(T, P)$, stability is determined by comparing $\ln \phi$:
   $$\Delta G = G_{\text{liq}} - G_{\text{vap}} = RT (\ln \phi_{\text{liq}} - \ln \phi_{\text{vap}})$$
   - If $\ln \phi_{\text{liq}} < \ln \phi_{\text{vap}} \iff \phi_{\text{liq}} < \phi_{\text{vap}}$: **Liquid phase is stable**.
   - If $\ln \phi_{\text{vap}} < \ln \phi_{\text{liq}} \iff \phi_{\text{vap}} < \phi_{\text{liq}}$: **Vapor phase is stable**.
   - If $\phi_{\text{liq}} = \phi_{\text{vap}}$: **Two-phase equilibrium ($P = P^{\text{sat}}(T)$)**.

### 6.3 Low-Pressure Regularization ($P \to 0$)
In the fugacity and residual property expressions, the term:
$$\xi(Z, B) = \frac{1}{2\sqrt{2} B} \ln \left( \frac{Z + (1 + \sqrt{2})B}{Z + (1 - \sqrt{2})B} \right)$$
takes the indeterminate form $\frac{0}{0}$ as $B \to 0$.

Rewrite the argument:
$$\frac{Z + (1 + \sqrt{2})B}{Z + (1 - \sqrt{2})B} = 1 + \frac{2\sqrt{2}B}{Z + (1 - \sqrt{2})B} = 1 + u$$
where $u = \frac{2\sqrt{2}B}{Z + (1 - \sqrt{2})B}$.
Using `std::log1p(u)`:
- For $u < 10^{-6}$, use Taylor expansion: $\ln(1+u) = u \left(1 - \frac{u}{2} + \frac{u^2}{3} - \frac{u^3}{4}\right)$, giving:
  $$\xi(Z, B) \approx \frac{1}{Z + (1-\sqrt{2})B} \left(1 - \frac{u}{2} + \frac{u^2}{3}\right)$$
- This avoids numerical cancellation and division by zero.

---

## 7. Numerical & Thermodynamic Pitfalls Identified

1. **Cardano Loss of Precision:** Using standard floating-point cube root on near-zero discriminants. *Mitigation:* Trigonometric branch for three real roots; clamped $\arccos$ arguments in $[-1, 1]$.
2. **Unphysical Intermediate Root Usage:** Evaluating properties at the spinodal saddle point. *Mitigation:* Explicit tagging of roots as `VAPOR`, `LIQUID`, `UNSTABLE_INTERMEDIATE`.
3. **Logarithm Domain Violations:** Evaluating $\ln(Z - B)$ when $Z \le B$. *Mitigation:* Pre-filter all candidate roots against $Z > B + \epsilon$.
4. **Unit Conversion Errors:** $P_c$ reported in bar or MPa while $R$ is in $\text{J}/(\text{mol}\cdot\text{K})$. *Mitigation:* Strict internal SI convention ($\text{Pa}, \text{K}, \text{m}^3/\text{mol}, \text{J}/\text{mol}$). User input/output conversion adapters.
5. **Phase Equilibrium Ambiguity:** Assuming largest root is always stable. *Mitigation:* Full fugacity/Gibbs energy comparison provided to the user.

---

## 8. Minimum Viable Feature Set (MVP) vs. Optional Extensions

### 8.1 MVP (Core Assignment Scope)
- Pure-component PR EOS (1976 formulation).
- Robust cubic solver (Cardano + Trigonometric).
- Forward PVT ($(T, P) \to Z, v$) with root classification and stable phase tagging.
- Inverse PVT ($(T, v) \to P$).
- Residual properties $H^{\text{res}}$ and $S^{\text{res}}$ with exact analytical $da/dT$.
- Fugacity coefficient $\phi$ and fugacity $f$ with small-$B$ regularization.
- Ideal gas heat capacity integration with standard reference state.
- Curated database of 8 common industrial gases.
- Interactive custom component input.
- Clean, menu-driven CLI interface.
- Complete automated unit test and validation suite.

### 8.2 Optional Future Enhancements (Post-MVP)
- Vapor pressure curve calculation via iterative flash solver ($f_{\text{liq}} = f_{\text{vap}}$).
- Peng–Robinson 1978 acentric factor correlation switch for $\omega > 0.49$.
- Export of PVT isotherms to CSV for external plotting.
