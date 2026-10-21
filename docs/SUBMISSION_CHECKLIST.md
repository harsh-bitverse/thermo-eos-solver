# CHL2002 Thermodynamics Course Project — Submission Checklist

**Project Title:** Prediction of Real-Gas PVT Behaviour, Thermodynamic Departure Functions, and Pure-Component Fugacity Using the Peng–Robinson Equation of State  
**Course:** CHL2002 Thermodynamics  
**Author:** Harsh Agarwal  
**Language & Standard:** C++17  
**Build System:** CMake 3.20+ (MinGW-W64 GCC 13.2)  
**Submission Due Date:** November 10, 2026  
**Repository Remote:** `https://github.com/harsh-bitverse/thermo-eos-solver.git`  
**Report Document:** `docs/FINAL_PROJECT_REPORT.pdf` (Strictly 4 Pages)  
**Source Archive:** `dist/thermo-eos-solver-source.zip`  

---

## 1. Compliance and Verification Checklist

| Item | Requirement / Criterion | Status | Verification Evidence / Reference |
| :--- | :--- | :---: | :--- |
| **1. Academic Deliverables** | 4-page technical report PDF (`FINAL_PROJECT_REPORT.pdf`) and clean computer program. | **VERIFIED** | PDF generated via ReportLab, verified strictly 4 pages. Verified with `pypdf`. |
| **2. C++17 Standard & Zero Dependencies** | Code strictly written in ISO C++17; relies exclusively on standard library headers. | **VERIFIED** | Built with `-std=c++17`; no external third-party dependencies required. |
| **3. Clean Compiler Warnings** | Zero compiler warnings under strict flags (`-Wall -Wextra -Wpedantic` / MSVC `/W4`). | **VERIFIED** | Clean build completed with 0 warnings under MinGW GCC 13.2.0. |
| **4. PR-1976 Formulation** | Canonical Peng–Robinson 1976 parameters ($a_c, b, \alpha(T), \kappa, da/dT$) correctly derived. | **VERIFIED** | Fully audited against Peng & Robinson (1976) Eqs. 11, 12, 13, 17, 18. |
| **5. Cubic Polynomial Solver** | Analytical Cardano / Vieta trigonometric solver; residual $|P(Z)| \le 10^{-12}$. | **VERIFIED** | Tested in `test_cubic_solver` (48 assertions), including critical point triple root and $\beta=1/3$ scaling. |
| **6. Mechanical Stability** | Analytical $(\partial P/\partial v)_T$ evaluation; classifies STABLE ($< 0$) vs UNSTABLE ($\ge 0$). | **VERIFIED** | Tested in `test_pvt_eos` (119 assertions) and `test_cli` (27 assertions). Discards unphysical roots ($Z \le B$). |
| **7. Residual Departure Functions** | Analytical $H^{\text{res}}$ and $S^{\text{res}}$ derived and implemented without empirical approximations. | **VERIFIED** | Verified against primary literature and independent differential finite differences. |
| **8. Pure Fugacity & ln(phi)** | Closed-form $\ln \phi$ and fugacity $f = \phi \cdot P$ implemented with small-$B$ Taylor regularization. | **VERIFIED** | Tested down to $P = 10^{-5}\text{ Pa}$; continuous to machine precision across threshold. |
| **9. Total Enthalpy & Entropy** | Reference state: ideal gas at $T_{\text{ref}} = 298.15\text{ K}$, $P_{\text{ref}} = 100\,000\text{ Pa}$ ($H=0, S=0$). | **VERIFIED** | NIST Shomate correlations integrated analytically; matches 1000-subinterval Simpson quadrature to $\le 10^{-6}$. |
| **10. Phase Preference Analysis** | Decides phase preference via Gibbs energy difference $\Delta G^{\text{res}} = R T \ln(\phi_L / \phi_V)$. | **VERIFIED** | Correctly resolves Propane at 320 K, 1 MPa ($\Delta G^{\text{res}} = +947.42\text{ J/mol} > 0 \implies$ Vapor preferred). |
| **11. Curated Gas Database** | Accurate properties ($T_c, P_c, \omega$) and Shomate intervals for 8 industrial gases. | **VERIFIED** | Methane, Ethane, Propane, n-Butane, Nitrogen, Oxygen, CO2, and Water Vapor verified against NIST SRD 69. |
| **12. Defensive Error Handling** | Rejects non-numeric, negative/zero inputs, and out-of-range temperatures without crashing. | **VERIFIED** | Water vapor steam reference limitation explicitly handled (`N/A`) with clear diagnostic message. |
| **13. Interactive CLI & Demos** | Easy-to-use menu with unit conversions and instant benchmark execution (`--demo all`). | **VERIFIED** | Built as `thermo_cli`; verified in `test_cli` and interactive execution. |
| **14. Automated Test Suite** | Comprehensive automated CTest suite with reproducible assertions. | **VERIFIED** | **373 / 373 assertions passing (100% pass rate)** across 5 CTest executables. |
| **15. True Differential Validation** | Avoids circular testing; verifies Maxwell and Gibbs–Helmholtz differential identities. | **VERIFIED** | Gibbs–Helmholtz: error $3.8 \times 10^{-10} \le 10^{-6}$; Pressure fugacity derivative: error $1.4 \times 10^{-9} \le 10^{-6}$. |
| **16. Clean Source Distribution** | Clean standalone source package excluding build directories, binaries, and caches. | **VERIFIED** | Packaged into `dist/thermo-eos-solver-source.zip` and verified building in a temporary directory. |
| **17. Mandatory Remote Sync** | Local branch `main` committed and synchronized with remote GitHub repository. | **VERIFIED** | Pushed to `origin/main` without force push; clean working tree. |

---

## 2. Directory Layout of Submission
```text
thermo-eos-solver/
├── CMakeLists.txt                      # Root CMake build configuration
├── README.md                           # Master user documentation and instructions
├── dist/
│   └── thermo-eos-solver-source.zip    # Clean source archive for distribution
├── docs/
│   ├── FINAL_PROJECT_REPORT.pdf        # Final 4-page academic technical report
│   ├── SUBMISSION_CHECKLIST.md         # This verification checklist
│   ├── PROJECT_SPECIFICATION.md        # Session 1 specification
│   ├── IMPLEMENTATION_ROADMAP.md       # Master project roadmap (Sessions 1-7 complete)
│   ├── SESSION_1_5_AUDIT.md            # Independent thermodynamic audit
│   ├── SESSION_3_RESULTS.md            # Session 3 results report
│   ├── SESSION_4_EQUATION_AUDIT.md     # Session 4 analytical derivations audit
│   ├── SESSION_4_RESULTS.md            # Session 4 results report
│   ├── SESSION_5_RESULTS.md            # Session 5 scientific validation report
│   ├── SESSION_6_RESULTS.md            # Session 6 interactive CLI report
│   └── SESSION_7_RESULTS.md            # Session 7 final audit & packaging report
├── include/
│   ├── cli/
│   │   └── CLIController.hpp           # Interactive CLI controller interface
│   └── thermo/
│       ├── Component.hpp               # Component & Shomate structures
│       ├── ComponentDatabase.hpp       # Pure component catalog interface
│       ├── Constants.hpp               # Physical constants & PR parameters
│       ├── CubicSolver.hpp             # Cardano / Vieta cubic solver interface
│       ├── PengRobinsonEOS.hpp         # PR-1976 equation of state engine
│       └── StateProperties.hpp         # PVTResult, PVTRoot, stability & phase tags
├── scripts/
│   └── generate_report_pdf.py          # ReportLab script for 4-page PDF report
├── src/
│   ├── cli/
│   │   └── CLIController.cpp           # Interactive CLI implementation
│   ├── cli_main.cpp                    # Binary entry point for thermo_cli
│   ├── Component.cpp                   # Component data & Shomate analytical integrals
│   ├── ComponentDatabase.cpp           # Curated database for 8 industrial gases
│   ├── CubicSolver.cpp                 # Analytical depressed cubic solver
│   └── PengRobinsonEOS.cpp             # Thermodynamics, departure functions & fugacity
└── tests/
    ├── CMakeLists.txt                  # Test target configuration
    ├── test_cubic_solver.cpp           # Cubic solver unit tests (48 assertions)
    ├── test_components.cpp             # Component & Simpson quadrature tests (44 assertions)
    ├── test_pvt_eos.cpp                # Forward/inverse PVT & reciprocity tests (119 assertions)
    ├── test_thermo_properties.cpp      # Departure functions & Maxwell tests (135 assertions)
    └── test_cli.cpp                    # Interactive CLI & recovery tests (27 assertions)
```

---

## 3. Fast Verification Commands for Evaluator

```bash
# 1. Clean configure and build
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build

# 2. Run the automated test suite (373 assertions, 5 targets)
ctest --test-dir build --output-on-failure

# 3. Run all three benchmark demonstration cases instantly
./build/thermo_cli --demo all

# 4. Launch interactive CLI
./build/thermo_cli
```
