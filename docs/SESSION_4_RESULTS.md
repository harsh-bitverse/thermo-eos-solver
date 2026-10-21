# CHL2002 Thermodynamics Project — Session 4: Residual Properties, Fugacity, Reference States, and Phase Equilibrium

**Document Version:** 1.0  
**Date:** October 21, 2026  
**Status:** Completed & Validated Module Report  

---

## 1. Executive Summary

In Session 4, the thermodynamic property and phase equilibrium layer of the CHL2002 pure-component thermodynamics engine was designed, implemented in C++17, and verified. Key achievements include:
1. **Residual Functions:** Implemented exact residual enthalpy $H^{\text{res}}(T, P, Z)$ and residual entropy $S^{\text{res}}(T, P, Z)$ utilizing analytical $\frac{da}{dT}$ and the temperature-dependent attractive parameter $a(T)$.
2. **Pure-Component Fugacity:** Implemented natural log fugacity coefficient $\ln \phi$, fugacity coefficient $\phi = \exp(\ln \phi)$, and fugacity $f = \phi P$ across all physically admissible roots ($Z > B$).
3. **Maxwell Identity Verification:** Proved and verified that the fundamental identity $S^{\text{res}} = \frac{H^{\text{res}}}{T} - R \ln \phi$ holds with absolute machine precision ($\le 1.0 \times 10^{-12}\text{ J}/(\text{mol}\cdot\text{K})$ error) across all tested states and roots.
4. **Low-Pressure Numerical Regularization:** Designed a dual-term Taylor regularization method for $u = \frac{2\sqrt{2}B}{Z + (1-\sqrt{2})B} < 10^{-5}$ and $\epsilon = Z - 1 - B < 10^{-5}$. This eliminates $0/0$ floating-point breakdown and catastrophic subtraction cancellation, ensuring smooth convergence to ideal-gas limits ($Z \to 1, \phi \to 1, f/P \to 1, H^{\text{res}} \to 0, S^{\text{res}} \to 0$) down to ultra-low pressures ($10^{-5}\text{ Pa}$).
5. **Reference-State Total Properties:** Integrated NIST Shomate ideal-gas integrals with the reference state $T_{\mathrm{ref}} = 298.15\text{ K}$, $P_{\mathrm{ref}} = 100\,000.0\text{ Pa}$. Verified that total $H(T_{\mathrm{ref}}, P_{\mathrm{ref}}) = 0.0\text{ J/mol}$ and $S(T_{\mathrm{ref}}, P_{\mathrm{ref}}) = 0.0\text{ J}/(\text{mol}\cdot\text{K})$ to within $10^{-10}$. Handled component temperature range limitations explicitly (e.g. Water Vapor with $T_{\mathrm{min}} = 500\text{ K} > T_{\mathrm{ref}}$ throws `std::out_of_range`).
6. **Thermodynamic Phase Equilibrium & Preference:** Classified multi-root states via $\Delta G_{\text{res}} = RT \ln(\phi_{\text{liq}} / \phi_{\text{vap}})$, strictly excluding mechanically unstable intermediate roots ($\partial P/\partial v > 0$). Verified subcooled liquid ($\phi_{\text{liq}} < \phi_{\text{vap}}$), superheated vapor ($\phi_{\text{vap}} < \phi_{\text{liq}}$), and liquid-vapor saturation coexistence ($P^{\text{sat}} \approx 1.603614\text{ MPa}$ for Propane at $320\text{ K}$).
7. **Comprehensive Testing:** 133 automated assertions in dedicated test target `test_thermo_properties` (100% passing), bringing the total test suite to 329 passing assertions with zero compiler warnings under `-Wall -Wextra -Wpedantic`.

---

## 2. Governing Equations & Numerical Conditioning

### 2.1 Thermodynamic Departure Functions
For a pure component described by the Peng–Robinson (1976) equation of state:
$$P = \frac{RT}{v - b} - \frac{a(T)}{v(v + b) + b(v - b)}$$
For any root satisfying $Z > B$ ($v > b$), define:
$$A = \frac{a(T) P}{(RT)^2}, \qquad B = \frac{b P}{RT}$$
$$L = \ln \left( \frac{Z + (1 + \sqrt{2})B}{Z + (1 - \sqrt{2})B} \right)$$

* **Residual Enthalpy:**
  $$H^{\text{res}}(T, P, Z) = RT(Z - 1) + \frac{T \frac{da}{dT} - a(T)}{2\sqrt{2}b} L \quad [\text{J}/\text{mol}]$$
* **Residual Entropy:**
  $$S^{\text{res}}(T, P, Z) = R \ln(Z - B) + \frac{\frac{da}{dT}}{2\sqrt{2}b} L \quad [\text{J}/(\text{mol}\cdot\text{K})]$$
* **Fugacity Coefficient:**
  $$\ln \phi = (Z - 1) - \ln(Z - B) - \frac{A}{2\sqrt{2}B} L, \qquad \phi = \exp(\ln \phi), \qquad f = \phi P \quad [\text{Pa}]$$

### 2.2 Numerical Strategy for Low Pressure ($B \to 0$)
Direct floating-point evaluation of $\frac{L}{2\sqrt{2}B}$ and $(Z - 1) - \ln(Z - B)$ can suffer catastrophic cancellation when $P \ll 1\text{ bar}$.

1. **Regularization of $L / (2\sqrt{2}B)$:**
   Let $u = \frac{2\sqrt{2}B}{Z + (1 - \sqrt{2})B}$. The logarithm argument is $1 + u$, so $L = \ln(1 + u)$.
   We define the shared dimensionless ratio:
   $$I_L = \frac{L}{2\sqrt{2}B} = \frac{1}{Z + (1 - \sqrt{2})B} \cdot \frac{\ln(1 + u)}{u}$$
   * When $|u| \ge 10^{-5}$: $I_L = \frac{\text{std::log1p}(u)}{2\sqrt{2}B}$, $L = \text{std::log1p}(u)$.
   * When $|u| < 10^{-5}$: 4-term Taylor expansion of $\frac{\ln(1+u)}{u}$:
     $$\frac{\ln(1+u)}{u} = 1 - \frac{u}{2} + \frac{u^2}{3} - \frac{u^3}{4} + \mathcal{O}(u^4)$$
     $$I_L = \frac{1}{Z + (1 - \sqrt{2})B} \left( 1 - \frac{u}{2} + \frac{u^2}{3} - \frac{u^3}{4} \right), \qquad L = 2\sqrt{2}B \cdot I_L$$
   * In $H^{\text{res}}$ and $S^{\text{res}}$, $\frac{L}{2\sqrt{2}b} = I_L \cdot \frac{P}{RT}$, which is unconditionally smooth and finite as $P \to 0$.

2. **Regularization of $(Z - 1) - \ln(Z - B)$:**
   Let $\epsilon = Z - 1 - B$. Then $Z - B = 1 + \epsilon$.
   $$(Z - 1) - \ln(1 + \epsilon) = (Z - 1) - \left( \epsilon - \frac{\epsilon^2}{2} + \frac{\epsilon^3}{3} - \frac{\epsilon^4}{4} \right) = B + \epsilon^2 \left( \frac{1}{2} - \frac{\epsilon}{3} + \frac{\epsilon^2}{4} \right)$$
   * When $|\epsilon| < 10^{-5}$: evaluate via $B + \epsilon^2 (0.5 - \epsilon/3 + \epsilon^2/4)$.
   * When $|\epsilon| \ge 10^{-5}$: evaluate via $(Z - 1.0) - \text{std::log}(Z - B)$.

### 2.3 Reference State Convention
* **Standard State:** $T_{\mathrm{ref}} = 298.15\text{ K}$, $P_{\mathrm{ref}} = 100\,000.0\text{ Pa}$ ($1.0\text{ bar}$).
* **Condition:** $H(T_{\mathrm{ref}}, P_{\mathrm{ref}}) = 0.0$, $S(T_{\mathrm{ref}}, P_{\mathrm{ref}}) = 0.0$.
* **Formulation:**
  $$H(T, P) = \Delta H^{\text{ig}}(T_{\mathrm{ref}} \to T) + H^{\text{res}}(T, P, Z) - H^{\text{res}}(T_{\mathrm{ref}}, P_{\mathrm{ref}}, Z_{\text{ref,vap}})$$
  $$S(T, P) = \Delta S^{\text{ig}}(T_{\mathrm{ref}} \to T) - R \ln\left(\frac{P}{P_{\mathrm{ref}}}\right) + S^{\text{res}}(T, P, Z) - S^{\text{res}}(T_{\mathrm{ref}}, P_{\mathrm{ref}}, Z_{\text{ref,vap}})$$
* **Note on Reference Root:** Departure terms at $(T_{\mathrm{ref}}, P_{\mathrm{ref}})$ are evaluated on the vapor-like root ($Z_{\text{ref,vap}} \approx 1$), which corresponds to the ideal-gas reference branch.

### 2.4 Thermodynamic Phase Preference
When multiple physical roots exist at $(T, P)$:
* The intermediate root is mechanically unstable ($\partial P/\partial v > 0$) and discarded from equilibrium comparisons.
* For the liquid-like root ($Z_{\min}$) and vapor-like root ($Z_{\max}$):
  $$\Delta G_{\text{res}} = G_{\text{liq}} - G_{\text{vap}} = RT \ln \left( \frac{\phi_{\text{liq}}}{\phi_{\text{vap}}} \right)$$
* **Classification Rules:**
  * If $\frac{|\phi_{\text{liq}} - \phi_{\text{vap}}|}{\phi_{\text{vap}}} \le 10^{-5} \implies$ `PhasePreference::EQUILIBRIUM` (saturation coexistence).
  * If $\phi_{\text{liq}} < \phi_{\text{vap}} \implies \Delta G_{\text{res}} < 0 \implies$ `PhasePreference::LIQUID_PREFERRED` (subcooled liquid).
  * If $\phi_{\text{vap}} < \phi_{\text{liq}} \implies \Delta G_{\text{res}} > 0 \implies$ `PhasePreference::VAPOR_PREFERRED` (superheated vapor).
  * If only 1 physical root exists $\implies$ `PhasePreference::SINGLE_PHASE` (no false phase coexistence).

---

## 3. Public API Additions

### 3.1 [`include/thermo/StateProperties.hpp`](file:///c:/Workspace/Projects/Thermo%20Software/include/thermo/StateProperties.hpp)
```cpp
namespace thermo {

enum class PhasePreference {
    SINGLE_PHASE,     ///< Exactly one physical root exists; no phase coexistence
    LIQUID_PREFERRED, ///< Multiple roots: phi_liq < phi_vap (G_liq < G_vap)
    VAPOR_PREFERRED,  ///< Multiple roots: phi_vap < phi_liq (G_vap < G_liq)
    EQUILIBRIUM       ///< Multiple roots: phi_liq == phi_vap within tolerance
};

struct PhaseEquilibriumResult {
    PhasePreference preference{PhasePreference::SINGLE_PHASE};
    double delta_G_res{0.0};  ///< Residual Gibbs difference: G_res,liq - G_res,vap = RT*ln(phi_liq / phi_vap) [J / mol]
    double phi_liq{0.0};      ///< Fugacity coefficient of liquid root
    double phi_vap{0.0};      ///< Fugacity coefficient of vapor root
    double phi_ratio{1.0};    ///< Fugacity coefficient ratio phi_liq / phi_vap
    bool is_equilibrium{false};///< True if |phi_liq - phi_vap| / phi_vap <= 1e-5
    std::string summary;      ///< Human-readable diagnostic description
};

struct PVTRoot {
    // ... PVT fields from Session 3 ...
    double H_res{0.0};           ///< Residual enthalpy in [J / mol]
    double S_res{0.0};           ///< Residual entropy in [J / (mol * K)]
    double ln_phi{0.0};          ///< Natural log of fugacity coefficient [-]
    double phi{1.0};             ///< Fugacity coefficient [-]
    double fugacity{0.0};        ///< Fugacity f = phi * P in [Pa]
    std::optional<double> H_total{std::nullopt}; ///< Total enthalpy relative to (T_ref, P_ref) [J / mol]
    std::optional<double> S_total{std::nullopt}; ///< Total entropy relative to (T_ref, P_ref) [J / (mol * K)]
};

struct PVTResult {
    // ... PVT fields from Session 3 ...
    PhaseEquilibriumResult equilibrium;

    [[nodiscard]] const PVTRoot* preferred_root() const noexcept;
};

} // namespace thermo
```

### 3.2 [`include/thermo/PengRobinsonEOS.hpp`](file:///c:/Workspace/Projects/Thermo%20Software/include/thermo/PengRobinsonEOS.hpp)
```cpp
namespace thermo {

class PengRobinsonEOS {
public:
    // ... Session 3 PVT API ...

    [[nodiscard]] double calculate_I_L(double Z, double B) const;
    [[nodiscard]] double calculate_L(double Z, double B) const;
    [[nodiscard]] double calculate_H_res(double T, double P, double Z) const;
    [[nodiscard]] double calculate_S_res(double T, double P, double Z) const;
    [[nodiscard]] double calculate_ln_phi(double T, double P, double Z) const;
    [[nodiscard]] double calculate_phi(double T, double P, double Z) const;
    [[nodiscard]] double calculate_fugacity(double T, double P, double Z) const;
    [[nodiscard]] double calculate_H_total(double T, double P, double Z) const;
    [[nodiscard]] double calculate_S_total(double T, double P, double Z) const;
    [[nodiscard]] PhaseEquilibriumResult evaluate_phase_preference(const PVTResult& pvt) const;
    [[nodiscard]] std::pair<double, double> reference_state_departures() const;
};

} // namespace thermo
```

---

## 4. Independent Benchmark Verification Table

All benchmark values were independently calculated and confirmed with zero discrepancies:

| Case | Substance | State ($T, P$) | Root | $Z$ | $H^{\text{res}}$ ($\text{J/mol}$) | $S^{\text{res}}$ ($\text{J}/(\text{mol}\cdot\text{K})$) | $\phi$ | $f$ ($\text{MPa}$) | Total $H$ ($\text{J/mol}$) | Total $S$ ($\text{J}/(\text{mol}\cdot\text{K})$) | Maxwell Identity Error |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **BM-01** | **Methane** | $300\text{ K}$, $5\text{ MPa}$ | Single | $0.901844$ | $-902.34$ | $-2.1446$ | $0.901385$ | $4.506926$ | $-818.18$ | $-34.408$ | $0.00 \times 10^{-16}$ |
| **BM-02** | **$\text{CO}_2$** | $300\text{ K}$, $4\text{ MPa}$ | Single | $0.752334$ | $-1990.99$ | $-4.7347$ | $0.795528$ | $3.182110$ | $-1881.04$ | $-35.083$ | $8.88 \times 10^{-16}$ |
| **BM-03** | **$\text{CO}_2$** | $300\text{ K}$, $7\text{ MPa}$ | Single | $0.198057$ | $-9374.37$ | $-27.5150$ | $0.638285$ | $4.467995$ | $-9264.43$ | $-62.516$ | $0.00 \times 10^{-16}$ |
| **BM-04** | **Propane** | $320\text{ K}$, $1\text{ MPa}$ | Liquid | $0.036453$ | $-14854.34$ | $-48.2035$ | $1.239280$ | $1.239280$ | $-12894.21$ | $-61.129$ | $0.00 \times 10^{-16}$ |
| | | | Intermediate | $0.091826$ | $-8210.60$ | $-28.3698$ | $1.385615$ | $1.385615$ | $-6250.47$ | $-41.296$ | $0.00 \times 10^{-16}$ |
| | | | Vapor | $0.850555$ | $-1131.04$ | $-2.3575$ | $0.868006$ | $0.868006$ | $+829.09$ | $-15.283$ | $4.44 \times 10^{-16}$ |
| **BM-05** | **Propane** | $320\text{ K}$, $1.8\text{ MPa}$ | Liquid | $0.063160$ | $-14545.92$ | $-47.2882$ | $0.708686$ | $1.275635$ | $-12585.79$ | $-65.101$ | $0.00 \times 10^{-16}$ |
| | | | Vapor | $0.704253$ | $-2355.85$ | $-5.9372$ | $0.764272$ | $1.375689$ | $-395.72$ | $-23.750$ | $0.00 \times 10^{-16}$ |
| **BM-06** | **Propane** | $320\text{ K}$, $P^{\text{sat}} \approx 1.603614\text{ MPa}$ | Liq / Vap | $0.057709$ / $0.736416$ | — | — | $\phi_{\text{liq}} = \phi_{\text{vap}} = 0.789884$ | $1.266670$ | — | — | Equilibrium ($\frac{\Delta \phi}{\phi} \approx 10^{-13}$) |

---

## 5. Test Suite Verification & Results

### CTest Summary
```
100% tests passed, 0 tests failed out of 4

Total Test time (real) = 0.03 sec
    Start 1: test_cubic_solver ................ Passed (41 assertions)
    Start 2: test_components .................. Passed (42 assertions)
    Start 3: test_pvt_eos ..................... Passed (113 assertions)
    Start 4: test_thermo_properties ........... Passed (133 assertions)
Total Assertion Count: 329 passed, 0 failed.
```

### Breakdown of the 133 New Assertions in `test_thermo_properties`:
* **Part A (Residual Identities, 26 checks):** Verified finite bounds, correct signs, and exact Maxwell identity $S^{\text{res}} = H^{\text{res}}/T - R \ln \phi$ ($\le 10^{-12}$) across Methane ($5\text{ MPa}$), $\text{CO}_2$ ($4\text{ MPa}$), and $\text{CO}_2$ ($7\text{ MPa}$).
* **Part B (Ideal-Gas Asymptotics, 36 checks):** Tested pressure sequences from $100\text{ kPa}$ down to $10^{-5}\text{ Pa}$. Verified $Z \to 1.0, \phi \to 1.0, H^{\text{res}} \to 0, S^{\text{res}} \to 0$. Verified continuity across the Taylor series threshold ($u = 10^{-5}$).
* **Part C (Fugacity & Three Roots, 22 checks):** Verified all 3 roots in Propane BM-04, mechanical instability of the intermediate root, and reciprocity $f = \phi P$.
* **Part D (Total H/S & Reference State, 25 checks):** Verified $H(T_{\mathrm{ref}}, P_{\mathrm{ref}}) = 0$ and $S(T_{\mathrm{ref}}, P_{\mathrm{ref}}) = 0$ for Methane and $\text{CO}_2$. Verified $\frac{\partial H}{\partial T} > 0$ and $\frac{\partial S}{\partial P} < 0$. Confirmed explicit `std::out_of_range` exception for Water Vapor ($T_{\mathrm{ref}} < 500\text{ K}$).
* **Part E (Phase Preference & Saturation, 17 checks):** Verified `VAPOR_PREFERRED` for Propane at $1.0\text{ MPa}$, `LIQUID_PREFERRED` at $1.8\text{ MPa}$, `EQUILIBRIUM` at $P^{\text{sat}} \approx 1.603614\text{ MPa}$, and `SINGLE_PHASE` for supercritical / dense single-root fluids.
* **Part F (Robustness & Domain Guardrails, 7 checks):** Rejection of $Z \le B$, $T \le 0$, $P \le 0$, NaN inputs, and out-of-range temperatures ($T > 1300\text{ K}$ for Methane).

---

## 6. Known Limitations

1. **Pure Components Only:** In strict compliance with project specifications, mixing rules and multicomponent flash algorithms are out of scope.
2. **Ideal-Gas Reference Range:** Total $H$ and $S$ rely on NIST Shomate ideal-gas integrals. If $T$ or $T_{\mathrm{ref}}$ lies outside $[T_{\mathrm{min}}, T_{\mathrm{max}}]$ (such as steam below $500\text{ K}$), total properties are withheld (`std::out_of_range`), while residual properties and fugacities remain fully operational.
3. **Model Nature of Cubic EOS:** Saturation pressure predictions are properties of the PR-1976 model and reflect analytical cubic behavior, distinct from empirical multiparameter Helmholtz equations of state.

---

## 7. Commands Used and Actual Results

1. **Build Core & Tests:**
   ```bash
   cmake --build build
   ```
   *Result:* Target `test_thermo_properties` compiled and linked cleanly. Zero compiler warnings under `-Wall -Wextra -Wpedantic`.

2. **Execute Dedicated Test:**
   ```bash
   .\build\tests\test_thermo_properties.exe
   ```
   *Result:* `SUMMARY: Passed 133 / 133 assertions.`

3. **Execute Full CTest Suite:**
   ```bash
   ctest --test-dir build --output-on-failure
   ```
   *Result:* `100% tests passed, 0 tests failed out of 4`.

---

## 8. Git Commit Record

* **Commit Hash:** `208198e`
* **Status:** Working tree clean, all 4 CTest targets passing (329/329 total assertions).

