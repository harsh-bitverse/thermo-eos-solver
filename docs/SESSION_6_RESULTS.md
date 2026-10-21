# CHL2002 Thermodynamics Project — Session 6: Interactive CLI, User Experience, and Demonstration Readiness

**Document Version:** 1.0  
**Date:** October 21, 2026  
**Status:** Completed & Synchronized Module Report  

---

## 1. Executive Summary

In Session 6, a modular, robust, interactive command-line interface (`thermo_cli`) was developed and integrated into the project. The interface enables students, evaluators, and engineers to easily demonstrate all core thermodynamic capabilities of the Peng–Robinson (PR-1976) software:

1. **Modular Architecture:** The CLI is implemented in `include/cli/CLIController.hpp` and `src/cli/CLIController.cpp`, strictly separated from the underlying thermodynamic calculations (`thermo_core`). The `CLIController` accepts configurable `std::istream` and `std::ostream` references, allowing full, automated, deterministic unit testing via standard memory streams.
2. **Interactive Menu Workflow:**
   - Option 1: Forward PVT behavior calculation and root stability classification.
   - Option 2: Thermodynamic properties (Residual Enthalpy/Entropy and Reference-State Total $H/S$).
   - Option 3: Pure-component fugacity and fugacity coefficient $\phi$.
   - Option 4: Comprehensive state analysis combining all thermodynamic properties into a unified report.
   - Option 5: Clean menu exit.
3. **Flexible Units with Transparent Conversions:**
   - Temperatures accepted in Kelvin [K] (SI default) or Celsius [°C] with explicit conversion ($T_{\text{K}} = T_{^\circ\text{C}} + 273.15$).
   - Pressures accepted in Pascal [Pa] (SI default), Bar [bar] ($1\text{ bar} = 100\,000\text{ Pa}$), or Megapascal [MPa] ($1\text{ MPa} = 1\,000\,000\text{ Pa}$).
4. **Input Guardrails & Exception Safety:**
   - Robust parsing rejects non-numeric input, negative/zero temperatures, and negative/zero pressures without throwing unhandled exceptions or crashing.
   - Immediate EOF or piped stream termination exits cleanly with code 0.
   - Temperature limits (such as the NIST Shomate $[500, 1700]\text{ K}$ range for Water Vapor where $T_{\text{ref}} = 298.15\text{ K}$ is out of range) are handled gracefully by displaying `N/A` with clear physical diagnostics rather than printing misleading zero values or crashing.
5. **Instant Benchmark Demonstration Mode:**
   - Added command-line demonstration flags (`--demo 1`, `--demo 2`, `--demo 3`, and `--demo all`) enabling an evaluator to reproduce the course benchmarks without manual keyboard entry.
6. **Automated Testing & Full Suite Health:**
   - Added `test_cli` (`tests/test_cli.cpp`) containing 27 targeted assertions covering input validation, recovery, units, Water Vapor limitations, 3-root phase preference, and direct binary invocation.
   - Total test suite expanded to **373 assertions across 5 CTest targets** with a **100% pass rate** and **zero compiler warnings** under `-Wall -Wextra -Wpedantic`.
7. **Git & Remote Synchronization:**
   - Synchronized and pushed to GitHub `origin/main` (`https://github.com/harsh-bitverse/thermo-eos-solver.git`).

---

## 2. Files Created and Modified

| Action | File Path | Purpose |
| :--- | :--- | :--- |
| **Created** | `include/cli/CLIController.hpp` | Header for `CLIController` defining interactive menus, prompt utilities, and formatted display routines. |
| **Created** | `src/cli/CLIController.cpp` | Implementation of stream-based interactive CLI, unit conversions, root formatting, and demonstration cases. |
| **Created** | `src/cli_main.cpp` | Main executable entry point for `thermo_cli`, processing command-line flags (`--help`, `--demo`). |
| **Created** | `tests/test_cli.cpp` | Comprehensive automated unit test suite for the CLI (27 assertions). |
| **Created** | `README.md` | Master project documentation, build instructions, CLI usage guide, and benchmark reference tables. |
| **Created** | `docs/SESSION_6_RESULTS.md` | Session 6 completion report and audit documentation. |
| **Modified** | `CMakeLists.txt` | Added `thermo_cli` executable target linked to `thermo_core`. |
| **Modified** | `tests/CMakeLists.txt` | Added `test_cli` CTest target. |
| **Modified** | `docs/IMPLEMENTATION_ROADMAP.md` | Marked Session 6 as completed. |

---

## 3. CLI Features and Usage

### Launching the Application
```bash
# Interactive mode:
./build/thermo_cli

# Automated demonstration mode:
./build/thermo_cli --demo all
```

### Demonstration Benchmark Outputs

#### Case 1: Methane at $300\text{ K}$, $5.0\text{ MPa}$ ($50\text{ bar}$)
```text
State Conditions:
  Temperature (T):     300.0000 K  (26.85 °C)  [Tr = 1.5743]
  Pressure (P):        5,000,000.00 Pa  (50.0000 bar = 5.0000 MPa)  [Pr = 1.0872]

[Root 1 - SINGLE-PHASE FLUID]
  Compressibility (Z): 0.901844  (Physical: Z > B = 0.053725)
  Molar Volume (v):    4.4990e-04 m^3/mol  (0.449901 L/mol)
  Isothermal (dP/dv)_T:-1.0153e+10 Pa/(m^3/mol)  [STABLE]
  Residual Enthalpy:   -902.34 J/mol
  Residual Entropy:    -2.1446 J/(mol*K)
  Fugacity Coeff (phi):0.901385  (ln(phi) = -0.103823)
  Fugacity (f):        4,506,925.70 Pa  (4.506926 MPa)
  Total Enthalpy (H):  -818.18 J/mol  [Ref: 298.15 K, 1.0 bar ideal gas]
  Total Entropy (S):   -34.4078 J/(mol*K) [Ref: 298.15 K, 1.0 bar ideal gas]

Phase Preference & Equilibrium Analysis:
  Thermodynamic Status: SINGLE_PHASE (Homogeneous single-phase fluid)
```

#### Case 2: Carbon Dioxide at $300\text{ K}$, $7.0\text{ MPa}$ ($70\text{ bar}$)
```text
State Conditions:
  Temperature (T):     300.0000 K  (26.85 °C)  [Tr = 0.9864]
  Pressure (P):        7,000,000.00 Pa  (70.0000 bar = 7.0000 MPa)  [Pr = 0.9489]

[Root 1 - SINGLE-PHASE FLUID]
  Compressibility (Z): 0.198057  (Physical: Z > B = 0.074837)
  Molar Volume (v):    7.0575e-05 m^3/mol  (0.070575 L/mol)
  Isothermal (dP/dv)_T:-8.8024e+10 Pa/(m^3/mol)  [STABLE]
  Residual Enthalpy:   -9374.37 J/mol
  Residual Entropy:    -27.5150 J/(mol*K)
  Fugacity Coeff (phi):0.638285  (ln(phi) = -0.448970)
  Fugacity (f):        4,467,995.09 Pa  (4.467995 MPa)
  Total Enthalpy (H):  -9264.43 J/mol  [Ref: 298.15 K, 1.0 bar ideal gas]
  Total Entropy (S):   -62.5164 J/(mol*K) [Ref: 298.15 K, 1.0 bar ideal gas]
```

#### Case 3: Propane at $320\text{ K}$, $1.0\text{ MPa}$ ($10\text{ bar}$)
```text
State Conditions:
  Temperature (T):     320.0000 K  (46.85 °C)  [Tr = 0.8653]
  Pressure (P):        1,000,000.00 Pa  (10.0000 bar = 1.0000 MPa)  [Pr = 0.2354]

[Root 1 - LIQUID-LIKE]
  Compressibility (Z): 0.036453  (Physical: Z > B = 0.021165)
  Molar Volume (v):    9.6988e-05 m^3/mol  (0.096988 L/mol)
  Isothermal (dP/dv)_T:-4.5722e+11 Pa/(m^3/mol)  [STABLE]
  Residual Enthalpy:   -14854.34 J/mol
  Residual Entropy:    -48.2035 J/(mol*K)
  Fugacity Coeff (phi):1.239280  (ln(phi) = 0.214531)
  Fugacity (f):        1,239,280.23 Pa  (1.239280 MPa)
  Total Enthalpy (H):  -12894.21 J/mol  [Ref: 298.15 K, 1.0 bar ideal gas]
  Total Entropy (S):   -61.1295 J/(mol*K) [Ref: 298.15 K, 1.0 bar ideal gas]

[Root 2 - INTERMEDIATE (UNPHYSICAL)]
  Compressibility (Z): 0.091826  (Physical: Z > B = 0.021165)
  Molar Volume (v):    2.4432e-04 m^3/mol  (0.244315 L/mol)
  Isothermal (dP/dv)_T: 1.8825e+10 Pa/(m^3/mol)  [UNSTABLE: (dP/dv)_T >= 0]
  Residual Enthalpy:   -8210.60 J/mol
  Residual Entropy:    -28.3698 J/(mol*K)
  Fugacity Coeff (phi):1.385615  (ln(phi) = 0.326144)
  Fugacity (f):        1,385,615.49 Pa  (1.385615 MPa)
  Notice:              Mechanically unstable root within spinodal region; unobservable.

[Root 3 - VAPOR-LIKE]
  Compressibility (Z): 0.850555  (Physical: Z > B = 0.021165)
  Molar Volume (v):    2.2630e-03 m^3/mol  (2.263012 L/mol)
  Isothermal (dP/dv)_T:-3.6879e+08 Pa/(m^3/mol)  [STABLE]
  Residual Enthalpy:   -1131.04 J/mol
  Residual Entropy:    -2.3575 J/(mol*K)
  Fugacity Coeff (phi):0.868006  (ln(phi) = -0.141556)
  Fugacity (f):        868,006.20 Pa  (0.868006 MPa)
  Total Enthalpy (H):   829.09 J/mol  [Ref: 298.15 K, 1.0 bar ideal gas]
  Total Entropy (S):   -15.2835 J/(mol*K) [Ref: 298.15 K, 1.0 bar ideal gas]

Phase Preference & Equilibrium Analysis:
  Thermodynamic Status: VAPOR_PREFERRED (Superheated Vapor)
  Preferred Root:       Vapor-like root (Z = 0.850555)
  Selection Criterion:  phi_vap (0.868006) < phi_liq (1.239280)
  Residual Gibbs Diff:  Delta G_res = G_liq - G_vap = +947.42 J/mol (> 0)
  Conclusion:           Vapor state has lower chemical potential and represents the stable phase.
  Thermodynamically Selected Root: Z = 0.850555 (VAPOR-LIKE)
```

---

## 4. Test Suite Audit & Assertion Breakdown

| Test Target Executable | Assertion Count | Scope & Focus |
| :--- | :---: | :--- |
| `test_cubic_solver` | **48** | Cardano/Vieta analytical cubic root finding, PR critical triple root, near-critical $\beta = 1/3$ scaling. |
| `test_components` | **44** | Component database lookups, Shomate integration, and 1000-subinterval Simpson quadrature validation. |
| `test_pvt_eos` | **119** | Forward/inverse PVT reciprocity, isothermal compressibility, critical inflection, $100\text{ MPa}$ supercritical extremes. |
| `test_thermo_properties` | **135** | Residual departure functions, pure fugacity, low-$P$ Taylor asymptotics, Maxwell and Gibbs–Helmholtz finite-difference checks. |
| `test_cli` | **27** | Stream-based CLI validation, input recovery, unit conversions, Water Vapor handling, 3-root phase preference, demonstration cases. |
| **Total Test Suite** | **373** | **100% Passed (0 Failures across 5 CTest targets)** |

---

## 5. Build Reproducibility & Compiler Cleanliness

Clean build executed from scratch:
```bash
cmake --build build --clean-first
```
- **Compiler:** MinGW-W64 GCC 13.2.0 (C++17)
- **Flags:** `-Wall -Wextra -Wpedantic`
- **Compiler Warnings:** **Zero warnings.**
- **CTest Execution:**
  ```text
  Test project C:/Workspace/Projects/Thermo Software/build
      Start 1: test_cubic_solver
  1/5 Test #1: test_cubic_solver ................   Passed    0.01 sec
      Start 2: test_components
  2/5 Test #2: test_components ..................   Passed    0.01 sec
      Start 3: test_pvt_eos
  3/5 Test #3: test_pvt_eos .....................   Passed    0.01 sec
      Start 4: test_thermo_properties
  4/5 Test #4: test_thermo_properties ...........   Passed    0.01 sec
      Start 5: test_cli
  5/5 Test #5: test_cli .........................   Passed    0.00 sec

  100% tests passed, 0 tests failed out of 5
  Total Test time (real) = 0.03 sec
  ```

---

## 6. Known Model Limitations & Physical Boundary Notes

1. **Pure Component Scope:** Only pure fluids are modeled; multicomponent mixtures and mixing rules are not supported.
2. **Water Vapor Standard State Limitation:** Gas-phase NIST Shomate parameters for Water Vapor are defined on $[500.0, 1700.0]\text{ K}$. Because the standard reference temperature $T_{\mathrm{ref}} = 298.15\text{ K}$ is outside this interval, total enthalpy and total entropy are explicitly reported as `N/A`. PVT, fugacity, and residual departure functions remain exact and fully functional.
3. **No Unphysical Extrapolations:** If an input temperature falls outside a component's documented Shomate interval, total property calculations are disallowed to prevent inaccurate polynomial extrapolation.

---

## 7. Git Synchronization

- **Target Branch:** `main`
- **Remote:** `origin` (`https://github.com/harsh-bitverse/thermo-eos-solver.git`)
- **Implementation Commit:** `14b60da`
- **Total Test Suite:** 373 assertions passing across 5 CTest targets (100% pass rate)
- **Status:** Synchronized cleanly with remote.
