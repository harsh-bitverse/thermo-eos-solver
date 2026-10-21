# CHL2002 Thermodynamics Project — Session 7: Final Audit, Technical Report & Submission Packaging

**Document Version:** 1.0  
**Date:** October 21, 2026  
**Status:** Completed, Verified, and Synchronized Final Report  

---

## 1. Executive Summary

Session 7 represents the completion of the CHL2002 Thermodynamics course project. All deliverables specified in the course mandate have been audited, finalized, verified, and packaged:

1. **Final Scientific & Numerical Audit:**
   - Every equation, parameter constant, unit conversion, and reference state was re-verified against canonical literature (Peng & Robinson 1976; NIST SRD 69; Poling et al. 2001).
   - All three benchmark demonstration cases (Methane at 300 K / 5 MPa, CO2 at 300 K / 7 MPa, and Propane at 320 K / 1 MPa) were directly reproduced from the compiled executable.
   - Reconciled the subcritical multi-root Propane Gibbs free energy difference to its exact IEEE 754 value ($\Delta G^{\text{res}} = +947.42\text{ J/mol} > 0$), confirming vapor phase preference.
2. **Academic Technical Report (`docs/FINAL_PROJECT_REPORT.pdf`):**
   - Produced a publication-quality technical report using ReportLab with custom running headers and footers.
   - **Page Budget:** Strictly **4 pages maximum** (verified as `Page 4 of 4` with zero overflow).
   - Formatted in professional academic style for a second-year chemical engineering thermodynamics course.
3. **Submission Verification Checklist (`docs/SUBMISSION_CHECKLIST.md`):**
   - Detailed 17-point compliance matrix proving adherence to every course requirement, physical limit, and architectural boundary.
4. **Clean Source Distribution Package (`dist/thermo-eos-solver-source.zip`):**
   - Packaged a clean standalone source archive containing 34 essential project files (sources, headers, CMake files, tests, documentation, and PDF report).
   - Excluded all build artifacts (`build/`, `.exe`, `.obj`, `.a`), git metadata, temporary files, and secrets.
   - Extracted into a clean temporary directory (`tmp_test_extract/`), built from scratch with CMake, and verified passing **100% of CTest targets** and demonstration commands.
5. **Full Test Suite & Compiler Health:**
   - **373 assertions passing across 5 CTest executables (100% pass rate, 0 failures)**.
   - Built under MinGW-W64 GCC 13.2.0 with `-Wall -Wextra -Wpedantic` with **zero warnings**.
6. **Mandatory GitHub Synchronization:**
   - Synchronized cleanly with remote `origin/main` (`https://github.com/harsh-bitverse/thermo-eos-solver.git`).

---

## 2. Final Benchmark Verification Table

Direct reproduction from `thermo_cli --demo all`:

| Benchmark Case | Substance & Conditions | Physical Roots | $Z$ | Molar Volume $v$ ($\text{m}^3/\text{mol}$) | Fugacity Coefficient $\phi$ | Residual Enthalpy $H^{\text{res}}$ ($\text{J/mol}$) | Residual Entropy $S^{\text{res}}$ ($\text{J}/(\text{mol}\cdot\text{K})$) | Mechanical Stability & Phase Verdict |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **Demo 1** | **Methane**<br>$300\text{ K}$, $5.0\text{ MPa}$ | 1 | $0.901844$ | $4.4990 \times 10^{-4}$ | $0.901385$<br>($f=4.507\text{ MPa}$) | $-902.34$ | $-2.1446$ | **STABLE**<br>`SINGLE_PHASE` (Supercritical gas) |
| **Demo 2** | **$\text{CO}_2$**<br>$300\text{ K}$, $7.0\text{ MPa}$ | 1 | $0.198057$ | $7.0575 \times 10^{-5}$ | $0.638285$<br>($f=4.468\text{ MPa}$) | $-9374.37$ | $-27.5150$ | **STABLE**<br>`SINGLE_PHASE` (Dense liquid-like fluid) |
| **Demo 3 (Root 1)** | **Propane**<br>$320\text{ K}$, $1.0\text{ MPa}$ | Root 1 (Liq) | $0.036453$ | $9.6988 \times 10^{-5}$ | $1.239280$ | $-14\,854.34$ | $-48.2035$ | **STABLE**<br>Metastable liquid branch |
| **Demo 3 (Root 2)** | **Propane**<br>$320\text{ K}$, $1.0\text{ MPa}$ | Root 2 (Mid) | $0.091826$ | $2.4432 \times 10^{-4}$ | $1.385615$ | $-8210.60$ | $-28.3698$ | **UNSTABLE** ($(\partial P/\partial v)_T > 0$)<br>Inside spinodal; unobservable |
| **Demo 3 (Root 3)** | **Propane**<br>$320\text{ K}$, $1.0\text{ MPa}$ | Root 3 (Vap) | $0.850555$ | $2.2630 \times 10^{-3}$ | $0.868006$ | $-1131.04$ | $-2.3575$ | **STABLE & PREFERRED**<br>$\Delta G^{\text{res}} = +947.42\text{ J/mol} > 0$ |

---

## 3. Automated Test Suite Final Summary

| Test Target | Executable | Assertion Count | Scope / Verification Focus |
| :--- | :--- | :---: | :--- |
| **1. Cubic Solver** | `test_cubic_solver` | **48** | Analytical Cardano/Vieta roots, critical triple root at $(T_c, P_c)$, and mean-field $\beta=1/3$ critical scaling. |
| **2. Components** | `test_components` | **44** | Database lookups, physical limits, and independent 1000-subinterval composite Simpson quadrature validation. |
| **3. PVT & EOS** | `test_pvt_eos` | **119** | Forward/inverse PVT reciprocity ($\le 10^{-10}$), mechanical stability, inflection at critical point, $100\text{ MPa}$ extremes. |
| **4. Thermodynamic Properties** | `test_thermo_properties` | **135** | Residual departure functions, pure fugacity, Taylor low-$P$ asymptotics, Gibbs–Helmholtz & Maxwell differential relations. |
| **5. CLI Controller** | `test_cli` | **27** | Stream-based menu navigation, input guardrails, unit conversions, Water Vapor handling, 3-root phase preference. |
| **Total Test Suite** | **All 5 Targets** | **373** | **100% Passed (373 / 373 assertions, 0 failed)** |

---

## 4. Technical Report Specifications

* **File Location:** `docs/FINAL_PROJECT_REPORT.pdf`
* **Page Count:** **4 pages exactly** (verified using `pypdf` and `pymupdf`).
* **Generation Engine:** Python `reportlab` (5.0.1) via `scripts/generate_report_pdf.py`.
* **Visual Verification:** Inspected with tool `view_file` — all typography, formulas, headers, footers, and tables fit without clipping.
* **Content Structure:**
  - Page 1: Title, metadata, introduction, physical scope, PR-1976 equations, analytical $da/dT$, cubic polynomial.
  - Page 2: Residual departure functions ($H^{\text{res}}, S^{\text{res}}$), fugacity ($\phi, f$), Shomate heat capacity integration, standard reference state ($298.15\text{ K}, 1.0\text{ bar}$), and low-$P$ Taylor series.
  - Page 3: Analytical Cardano/Vieta solver, mechanical stability criterion ($(\partial P/\partial v)_T < 0$), phase preference criterion, C++17 modular architecture table, and compilation commands.
  - Page 4: Audited benchmark verification table, in-depth Propane multi-root analysis, test suite and Maxwell cross-check summary, assumptions, model limitations, conclusions, and primary literature references.

---

## 5. Source Distribution Archive Verification

* **Archive Location:** `dist/thermo-eos-solver-source.zip`
* **Total Packaged Files:** 34 files (size: $\sim 128\text{ KB}$).
* **Build Artifact Filtering:** Zero build artifacts, zero binaries, zero caches.
* **Clean Build Verification:** Extracted into a temporary directory, configured with CMake, compiled with GCC, and executed CTest:
  ```text
  100% tests passed, 0 tests failed out of 5
  thermo_cli --demo 1 passed in extracted source package!
  ```

---

## 6. Build Environment and Compiler Status

* **Operating System:** Windows (MinGW-W64 GCC 13.2.0)
* **Standard:** C++17 (`-std=c++17`)
* **Flags:** `-Wall -Wextra -Wpedantic`
* **Warnings:** **Zero compiler warnings.**

---

## 7. Submission Package File Inventory

* `CMakeLists.txt`
* `README.md`
* `dist/thermo-eos-solver-source.zip`
* `docs/FINAL_PROJECT_REPORT.pdf`
* `docs/SUBMISSION_CHECKLIST.md`
* `docs/IMPLEMENTATION_ROADMAP.md`
* `docs/PROJECT_SPECIFICATION.md`
* `docs/SESSION_1_5_AUDIT.md`
* `docs/SESSION_3_RESULTS.md`
* `docs/SESSION_4_EQUATION_AUDIT.md`
* `docs/SESSION_4_RESULTS.md`
* `docs/SESSION_5_RESULTS.md`
* `docs/SESSION_6_RESULTS.md`
* `docs/SESSION_7_RESULTS.md`
* `include/cli/CLIController.hpp`
* `include/thermo/Component.hpp`
* `include/thermo/ComponentDatabase.hpp`
* `include/thermo/Constants.hpp`
* `include/thermo/CubicSolver.hpp`
* `include/thermo/PengRobinsonEOS.hpp`
* `include/thermo/StateProperties.hpp`
* `scripts/generate_report_pdf.py`
* `scripts/package_source.py`
* `src/cli_main.cpp`
* `src/cli/CLIController.cpp`
* `src/Component.cpp`
* `src/ComponentDatabase.cpp`
* `src/CubicSolver.cpp`
* `src/PengRobinsonEOS.cpp`
* `tests/CMakeLists.txt`
* `tests/test_cubic_solver.cpp`
* `tests/test_components.cpp`
* `tests/test_pvt_eos.cpp`
* `tests/test_thermo_properties.cpp`
* `tests/test_cli.cpp`

---

## 8. Git Commit and Remote Synchronization

* **Repository:** `https://github.com/harsh-bitverse/thermo-eos-solver.git`
* **Branch:** `main` tracking `origin/main`
* **Commit:** `8d271d3`
* **Status:** Synchronized cleanly with remote.

