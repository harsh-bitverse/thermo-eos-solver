# CHL2002 Thermodynamics Project — Session 8: Scientific Audit & Verification Record

**Document Version:** 1.0  
**Date:** November 2026 (Audited October 21, 2026)  
**Author:** Harshith N  
**Course:** CHL2002 Thermodynamics  
**Status:** Complete & Verified  

---

## 1. Executive Summary of Audit

Session 8 performed a targeted scientific and numerical audit of the pure-component Peng–Robinson (PR-1976) software implementation, the 4-page academic technical report (`docs/FINAL_PROJECT_REPORT.pdf`), user documentation (`README.md`), and submission checklist (`docs/SUBMISSION_CHECKLIST.md`).

All factual claims were cross-referenced against the actual C++17 source code, unit tests, and primary scientific references. Discrepancies between historical session reports and the actual codebase were identified and reconciled.

---

## 2. Files Inspected

1. **Core Thermodynamics Library (`thermo_core`):**
   - `include/thermo/Constants.hpp` (Universal constants, analytical PR parameters, numerical tolerances)
   - `include/thermo/Component.hpp` and `src/Component.cpp` (PureComponent structure, validation, Shomate integration)
   - `include/thermo/ComponentDatabase.hpp` and `src/ComponentDatabase.cpp` (Curated 8-fluid catalog)
   - `include/thermo/CubicSolver.hpp` and `src/CubicSolver.cpp` (Cardano and Vieta trigonometric polynomial solver)
   - `include/thermo/StateProperties.hpp` (PVTRoot, PVTResult, PhasePreference, stability tags)
   - `include/thermo/PengRobinsonEOS.hpp` and `src/PengRobinsonEOS.cpp` (PR engine, departure functions, fugacity, phase preference)
2. **User Interface Module:**
   - `include/cli/CLIController.hpp` and `src/cli/CLIController.cpp` (Stream-decoupled interactive menu & demo runners)
   - `src/cli_main.cpp` (Binary entry point and CLI argument parser)
3. **Automated Test Targets:**
   - `tests/test_cubic_solver.cpp` (48 assertions)
   - `tests/test_components.cpp` (44 assertions)
   - `tests/test_pvt_eos.cpp` (119 assertions)
   - `tests/test_thermo_properties.cpp` (135 assertions)
   - `tests/test_cli.cpp` (27 assertions)
4. **Documentation & Build Configuration:**
   - `CMakeLists.txt` and `tests/CMakeLists.txt`
   - `README.md`
   - `scripts/generate_report_pdf.py` and `docs/FINAL_PROJECT_REPORT.pdf`
   - `docs/SUBMISSION_CHECKLIST.md`

---

## 3. Verified Report Claims Against Source Code & Tests

| Claim in Report | Source Code / Test Verification | Verification Status |
| :--- | :--- | :---: |
| **Universal Gas Constant** | `constants::R = 8.31446261815324 J/(mol*K)` in `Constants.hpp`. | **VERIFIED** (CODATA 2018 exact defining constant) |
| **PR Parameters Precision** | `OMEGA_A = 0.4572355289213822`, `OMEGA_B = 0.07779607390388846` in `Constants.hpp`. | **VERIFIED** (Derived from cubic inflection constraints rather than 5-decimal rounded textbook values) |
| **Analytical $da/dT$** | `PengRobinsonEOS::da_dT(double T)` in `PengRobinsonEOS.cpp:92`. | **VERIFIED** (Matches central difference to relative error $\le 10^{-7}$ in `test_pvt_eos`) |
| **Cubic Solver Residual** | `CubicSolver::solve()` guarantees $|P(Z)| \le 10^{-12}$ for regular roots. | **VERIFIED** (Tested across 48 cubic test cases in `test_cubic_solver`) |
| **Root Filtering $Z > B$** | `PengRobinsonEOS.cpp:196` strictly filters $Z > B + 10^{-12}$. | **VERIFIED** (Discards unphysical sub-covolume roots) |
| **Mechanical Stability** | `PengRobinsonEOS::calculate_dP_dv(T, v)` evaluated analytically. | **VERIFIED** (STABLE if $< 0$, UNSTABLE if $\ge 0$; tested in `test_pvt_eos`) |
| **Taylor Regularization** | `calculate_I_L()` uses $u = 2\sqrt{2}B / [Z + (1-\sqrt{2})B]$ with 4-term expansion for $|u| < 10^{-5}$. | **VERIFIED** (Continuous to machine precision down to $P = 10^{-5}\text{ Pa}$) |
| **Reference State** | $T_{\text{ref}} = 298.15\text{ K}$, $P_{\text{ref}} = 100\,000\text{ Pa}$ ($1.0\text{ bar}$). | **VERIFIED** (Total $H = 0, S = 0$ at reference state in `test_thermo_properties`) |
| **Maxwell & Gibbs–Helmholtz** | Differential cross-checks via central finite differences. | **VERIFIED** (Gibbs–Helmholtz error $3.8 \times 10^{-10} \le 10^{-6}$; $\partial \ln \phi/\partial P$ error $1.4 \times 10^{-9} \le 10^{-6}$) |
| **Water Vapor Limitation** | Throws `std::out_of_range` for total $H/S$ because $298.15\text{ K} < 500.0\text{ K}$. | **VERIFIED** (Tested in `test_thermo_properties` and `test_cli`; reported as `N/A`) |
| **Propane Multi-Root State** | At 320 K, 1 MPa: 3 real roots; Root 2 unstable ($dP/dv > 0$); Root 3 vapor preferred ($\Delta G^{\text{res}} = +947.42\text{ J/mol}$). | **VERIFIED** (Reproduced exactly in binary demo and `test_thermo_properties`) |

---

## 4. Discrepancies Found and Resolutions

1. **Author Name Discrepancy:**
   - *Finding:* Early draft documents used the placeholder name “Harsh Agarwal”.
   - *Resolution:* The student explicitly confirmed the author name as **Harshith N**. Updated in `scripts/generate_report_pdf.py`, `docs/FINAL_PROJECT_REPORT.pdf`, and `docs/SUBMISSION_CHECKLIST.md`.
2. **Propane Phase Envelope Terminology (Page 4):**
   - *Finding:* The report originally described Propane at 320 K, 1.0 MPa as lying in the “subcritical two-phase region”.
   - *Correction:* Because 1.0 MPa is below the saturation pressure ($P^{\text{sat}} \approx 1.6036\text{ MPa}$), liquid and vapor do not coexist in equilibrium. The state is in the superheated vapor regime. Replaced with:
     `"subcritical three-root region of the Peng–Robinson equation"`,
     explaining that the existence of three physical roots reflects the mathematical presence of spinodal branches and does not alone prove two-phase coexistence.
3. **Propane $\Delta G^{\text{res}}$ Numerical Value:**
   - *Finding:* In some earlier notes, $\Delta G^{\text{res}}$ was quoted as $+947.01\text{ J/mol}$ from lower-precision manual calculations, whereas the compiled program computes $+947.42\text{ J/mol}$.
   - *Resolution:* Updated `CLIController.cpp` line 738 and all documentation to record the exact IEEE 754 value $+947.42\text{ J/mol}$.
4. **Literature-Verification Claim (Page 4):**
   - *Finding:* The report stated calculations were “independently verified against literature benchmarks and internal mathematical consistency checks”. The automated tests test PR exact equations rather than experimental data tables.
   - *Correction:* Replaced with the precise, technically defensible statement:
     `"The implementation was evaluated using direct benchmark demonstrations and independent mathematical consistency checks."`
5. **Taylor Regularization Formula in Report (Page 2):**
   - *Finding:* The report originally stated a series in powers of $(B/Z)$, whereas the C++ code implements $\ln(1+u)/u$ where $u = 2\sqrt{2}B / [Z + (1-\sqrt{2})B]$.
   - *Correction:* Updated equation in `scripts/generate_report_pdf.py` and `FINAL_PROJECT_REPORT.pdf` to match the exact formula in `PengRobinsonEOS.cpp:306-313`.
6. **Phase Coexistence vs. Saturation Pressure Solving (Page 3):**
   - *Finding:* The report stated that if $|\phi_{\text{liq}} - \phi_{\text{vap}}| / \phi_{\text{vap}} \le 10^{-5}$, equilibrium exists, without clarifying that the program evaluates state properties at specified inputs rather than searching for unknown saturation pressures.
   - *Correction:* Explicitly qualified that the program identifies approximate coexistence at user-provided states (e.g. BM-06 at 1.6036 MPa) rather than running an iterative VLE bubble/dew-point algorithm from scratch.

---

## 5. Build and Test Verification Evidence

```text
Compilation Command: cmake -B build -S . -DCMAKE_BUILD_TYPE=Release && cmake --build build
Compiler: MinGW-W64 GCC 13.2.0 (C++17)
Flags: -Wall -Wextra -Wpedantic
Warnings: Zero warnings (0 warnings across all source units)

Test Suite Execution: ctest --test-dir build --output-on-failure
Output:
    Start 1: test_cubic_solver
1/5 Test #1: test_cubic_solver ................   Passed    0.01 sec  (48 assertions)
    Start 2: test_components
2/5 Test #2: test_components ..................   Passed    0.01 sec  (44 assertions)
    Start 3: test_pvt_eos
3/5 Test #3: test_pvt_eos .....................   Passed    0.01 sec  (119 assertions)
    Start 4: test_thermo_properties
4/5 Test #4: test_thermo_properties ...........   Passed    0.01 sec  (135 assertions)
    Start 5: test_cli
5/5 Test #5: test_cli .........................   Passed    0.01 sec  (27 assertions)

100% tests passed, 0 tests failed out of 5
Total Test time (real) = 0.04 sec
Total Test Assertions: 373 passed, 0 failed
```

---

## 6. Key Viva Concepts the Student Must Be Prepared to Explain

1. **Why three roots exist in the cubic equation of state:**  
   Subcritical isotherms exhibit a van der Waals "S-curve" loop between the liquid and vapor spinodals. The cubic polynomial in $Z$ has three real roots ($\Delta < 0$). The smallest root is liquid-like, the largest root is vapor-like, and the intermediate root lies along the spinodal slope where $(\partial P/\partial v)_T > 0$.
2. **Why the intermediate root is physically unobservable:**  
   Mechanical stability requires that compressing a fluid increases its pressure: $(\partial P/\partial v)_T < 0$. If $(\partial P/\partial v)_T > 0$, an infinitesimal density fluctuation spontaneously amplifies (spinodal decomposition), making the state unobservable in a static experiment.
3. **How phase preference is determined:**  
   At fixed $T$ and $P$, the stable phase is the one that minimizes Gibbs free energy $G$. Because $G_i = G^{\text{ig}} + R T \ln \phi_i$, comparing Gibbs energies is mathematically equivalent to comparing fugacity coefficients:
   $$\Delta G^{\text{res}} = G_{\text{liq}} - G_{\text{vap}} = R T \ln\left(\frac{\phi_{\text{liq}}}{\phi_{\text{vap}}}\right)$$
   The phase with the lower fugacity coefficient has lower chemical potential and is preferred.
4. **Difference between Phase Preference and Saturation Calculation:**  
   The program evaluates the properties of liquid and vapor branches at the given input pressure. It does not iterate pressure using a numerical solver to find where $\phi_{\text{liq}} = \phi_{\text{vap}}$ (saturation pressure). When given the saturation pressure (e.g. 1.6036 MPa for Propane at 320 K), it confirms $\phi_{\text{liq}} \approx \phi_{\text{vap}}$.
5. **The Steam (Water Vapor) Shomate Limitation:**  
   Gas-phase NIST Shomate heat capacity coefficients for Water Vapor are valid only on $[500.0, 1700.0]\text{ K}$. Because the standard reference temperature $T_{\text{ref}} = 298.15\text{ K}$ lies below $500\text{ K}$, integrating the correlation down to $T_{\text{ref}}$ would constitute unphysical polynomial extrapolation. Hence, total $H$ and $S$ relative to 298.15 K are reported as `N/A`. PVT, fugacity, and residual departure properties remain exact.

---

## 7. Viva Voce Preparation Guide Created

To ensure comprehensive preparation for defense before course evaluators, a dedicated 16-topic study guide was created:
- **Path:** [`docs/VIVA_PREPARATION.md`](file:///c:/Workspace/Projects/Thermo%20Software/docs/VIVA_PREPARATION.md)
- **Contents:**
  1. *Core Thermodynamics Questions & Answers (Q1–Q10):* Non-ideal gas physics, physical role of $a(T)$ and $b$, Pitzer acentric factor, derivation of PR cubic polynomial, departure functions, fugacity and chemical potential equilibrium, mechanical stability $(\partial P/\partial v)_T < 0$, Gibbs energy minimization ($\Delta G^{\text{res}} = R T \ln(\phi_L / \phi_V)$), state evaluation vs. saturation pressure solving, and steam reference limitation.
  2. *Numerical and Implementation Questions (Q11–Q16):* Tschirnhaus depression, Cardano formula vs. Vieta trigonometric method for $\Delta \le 0$, root filtering ($Z > B \iff v > b$), exact analytical $da/dT$ vs. finite difference noise, Taylor series expansion for $u = 2\sqrt{2}B / [Z + (1-\sqrt{2})B]$ below $|u| < 10^{-5}$, and precise test verification scope.
  3. *Code Walkthrough & Trace:* Equation-to-file-and-function architecture mapping table, followed by an end-to-end trace from binary launch to multi-root phase resolution.

