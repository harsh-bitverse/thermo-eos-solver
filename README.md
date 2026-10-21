# CHL2002 Thermodynamics Software: Peng–Robinson EOS Solver

A modern, scientifically validated C++17 engineering application for pure-component real-gas thermodynamic property calculations using the canonical Peng–Robinson (PR-1976) cubic equation of state.

Developed for the **CHL2002 Thermodynamics** course.

---

## 1. Project Overview & Scope

This software predicts the pressure–volume–temperature (PVT) behavior and thermodynamic departure properties of pure substances. The architecture strictly separates core thermodynamic calculations, numerical cubic polynomial solvers, curated property databases, and an interactive command-line interface (CLI).

### Core Capabilities
1. **PVT State Behavior:**
   - Evaluates forward PVT at $(T, P)$ by solving the dimensionless cubic polynomial in compressibility factor $Z$.
   - Discards unphysical roots ($Z \le B$).
   - Evaluates the analytical isothermal derivative $(\partial P/\partial v)_T$ to classify mechanical stability (STABLE vs. UNSTABLE).
   - Computes inverse PVT: $P(T, v)$ directly from the Peng–Robinson equation of state.
2. **Residual Departure Functions:**
   - Residual Enthalpy $H^{\text{res}}(T, P, Z) = H(T, P) - H^{\text{ig}}(T, P)$ in $\text{J/mol}$.
   - Residual Entropy $S^{\text{res}}(T, P, Z) = S(T, P) - S^{\text{ig}}(T, P)$ in $\text{J}/(\text{mol}\cdot\text{K})$.
   - Regularized small-$B$ Taylor expansions prevent numerical breakdown at ultra-low pressures ($P < 10^{-5}\text{ Pa}$).
3. **Pure-Component Fugacity:**
   - Natural log of fugacity coefficient $\ln \phi$ and fugacity coefficient $\phi = \exp(\ln \phi)$.
   - Pure-component fugacity $f = \phi \cdot P$ in $\text{Pa}$ and $\text{MPa}$.
4. **Ideal-Gas Heat Capacity & Total Properties:**
   - Authentic 5-parameter NIST Shomate correlations:
     $$C_p^\circ(t) = A + B t + C t^2 + D t^3 + \frac{E}{t^2} \quad (t = T / 1000)$$
   - Analytical integration of $\Delta H^{\text{ig}}(T_{\text{ref}} \to T)$ and $\Delta S^{\text{ig}}(T_{\text{ref}} \to T)$.
   - Standard reference-state total enthalpy $H(T, P)$ and total entropy $S(T, P)$ referenced to:
     $$T_{\text{ref}} = 298.15\text{ K} \; (25\,^\circ\text{C}), \quad P_{\text{ref}} = 100\,000\text{ Pa} \; (1.0\text{ bar})$$
     where $H_{\text{total}}(T_{\text{ref}}, P_{\text{ref}}) \equiv 0.0\text{ J/mol}$ and $S_{\text{total}}(T_{\text{ref}}, P_{\text{ref}}) \equiv 0.0\text{ J}/(\text{mol}\cdot\text{K})$.
5. **Phase Preference & Equilibrium Analysis:**
   - In three-root regions, evaluates Gibbs energy difference $\Delta G^{\text{res}} = R T \ln(\phi_L / \phi_V)$.
   - Determines thermodynamically stable phase (`LIQUID_PREFERRED`, `VAPOR_PREFERRED`, or `EQUILIBRIUM`).
   - Automatically identifies and isolates mechanically unstable intermediate roots within the spinodal zone.

---

## 2. Supported Pure Components

The database provides curated physical parameters ($T_c, P_c, \omega$, molar mass) and authentic NIST Shomate parameters for 8 industrial fluids:

| Index | Substance | Formula | Molar Mass ($\text{kg/mol}$) | $T_c$ ($\text{K}$) | $P_c$ ($\text{MPa}$) | Acentric Factor $\omega$ | Shomate Validity Range |
| :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **1** | Methane | $\text{CH}_4$ | 0.016043 | 190.56 | 4.599 | 0.0115 | $[298.0, 1300.0]\text{ K}$ |
| **2** | Ethane | $\text{C}_2\text{H}_6$ | 0.030070 | 305.32 | 4.872 | 0.0995 | $[298.0, 1200.0]\text{ K}$ |
| **3** | Propane | $\text{C}_3\text{H}_8$ | 0.044096 | 369.83 | 4.248 | 0.1523 | $[298.0, 1200.0]\text{ K}$ |
| **4** | n-Butane | $\text{C}_4\text{H}_{10}$ | 0.058122 | 425.12 | 3.796 | 0.2002 | $[298.0, 1200.0]\text{ K}$ |
| **5** | Nitrogen | $\text{N}_2$ | 0.028013 | 126.20 | 3.398 | 0.0372 | $[100.0, 500.0]\text{ K}$ |
| **6** | Oxygen | $\text{O}_2$ | 0.031999 | 154.58 | 5.043 | 0.0222 | $[100.0, 700.0]\text{ K}$ |
| **7** | Carbon Dioxide | $\text{CO}_2$ | 0.044010 | 304.13 | 7.377 | 0.2239 | $[298.0, 1200.0]\text{ K}$ |
| **8** | Water Vapor | $\text{H}_2\text{O}$ | 0.018015 | 647.10 | 22.064 | 0.3443 | $[500.0, 1700.0]\text{ K}$ |

---

## 3. Build & Requirements

### Prerequisites
- C++17 compliant compiler (GCC/MinGW $\ge 9.0$, Clang $\ge 10.0$, or MSVC $\ge 2019$)
- CMake $\ge 3.20$
- Standard C++ library only (no external dependencies)

### Build Instructions
```bash
# Configure build directory
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release

# Compile all targets (core library, CLI, and test suite)
cmake --build build

# Run automated test suite (373 assertions across 5 CTest targets)
ctest --test-dir build --output-on-failure
```

---

## 4. CLI Usage & Features

Launch the interactive interface:
```bash
./build/thermo_cli
```

### Main Menu Options
```text
MAIN MENU:
  [1] Calculate PVT properties for a pure component
  [2] Calculate thermodynamic properties H and S (Residual & Total)
  [3] Calculate pure-component fugacity and fugacity coefficient
  [4] Display all available properties for one state (Combined)
  [5] Exit
```

### Flexible Unit Support
- **Temperature:** Accepts values in **Kelvin [K]** or **Celsius [°C]** (with explicit conversion displayed).
- **Pressure:** Accepts values in **Pascal [Pa]**, **Bar [bar]**, or **Megapascal [MPa]**.

### Robust Input Guardrails
- Automatically validates component numbers, names, or formulas (case-insensitive).
- Rejects non-numeric text, negative or zero temperatures, and negative or zero pressures without crashing.
- Handles unexpected EOF cleanly without infinite loops.
- Gracefully explains temperature limitations when outside NIST Shomate intervals.

---

## 5. Benchmark Demonstration Cases

The software includes three reproducible benchmark cases from the CHL2002 course specification. You can run them instantly via command-line flags or through the interactive menu:

### Command-Line Execution
```bash
# Run all three demonstrations
./build/thermo_cli --demo all

# Or run individual cases
./build/thermo_cli --demo 1    # Methane (300 K, 5 MPa)
./build/thermo_cli --demo 2    # CO2 (300 K, 7 MPa)
./build/thermo_cli --demo 3    # Propane (320 K, 1 MPa)
```

### Benchmark Reference Results

#### Demo 1: Methane at $300\text{ K}$, $5.0\text{ MPa}$ ($50\text{ bar}$) — Supercritical Single-Phase Gas
- **State Classification:** Single-phase fluid ($\Delta > 0$, 1 root)
- **Compressibility Factor $Z$:** $0.901844$
- **Molar Volume $v$:** $4.4990 \times 10^{-4}\text{ m}^3/\text{mol}$ ($0.4499\text{ L/mol}$)
- **Fugacity Coefficient $\phi$:** $0.901385$ ($\ln \phi = -0.103823$)
- **Fugacity $f$:** $4.506926\text{ MPa}$ ($4\,506\,926\text{ Pa}$)
- **Residual Enthalpy $H^{\text{res}}$:** $-902.34\text{ J/mol}$
- **Residual Entropy $S^{\text{res}}$:** $-2.1446\text{ J}/(\text{mol}\cdot\text{K})$
- **Total Enthalpy $H$:** $-818.18\text{ J/mol}$
- **Total Entropy $S$:** $-34.4078\text{ J}/(\text{mol}\cdot\text{K})$
- **Stability:** Mechanically stable ($(\partial P/\partial v)_T = -1.0153 \times 10^{10}\text{ Pa}/(\text{m}^3/\text{mol}) < 0$)

#### Demo 2: Carbon Dioxide at $300\text{ K}$, $7.0\text{ MPa}$ ($70\text{ bar}$) — Dense Compressed Fluid
- **State Classification:** Single-phase compressed liquid-like fluid near critical point ($T_c = 304.13\text{ K}$)
- **Compressibility Factor $Z$:** $0.198057$
- **Molar Volume $v$:** $7.0575 \times 10^{-5}\text{ m}^3/\text{mol}$ ($0.07058\text{ L/mol}$)
- **Fugacity Coefficient $\phi$:** $0.638285$ ($\ln \phi = -0.448970$)
- **Fugacity $f$:** $4.467995\text{ MPa}$ ($4\,467\,995\text{ Pa}$)
- **Residual Enthalpy $H^{\text{res}}$:** $-9374.37\text{ J/mol}$
- **Residual Entropy $S^{\text{res}}$:** $-27.5150\text{ J}/(\text{mol}\cdot\text{K})$
- **Total Enthalpy $H$:** $-9264.43\text{ J/mol}$
- **Total Entropy $S$:** $-62.5164\text{ J}/(\text{mol}\cdot\text{K})$
- **Stability:** Mechanically stable ($(\partial P/\partial v)_T = -8.8024 \times 10^{10}\text{ Pa}/(\text{m}^3/\text{mol}) < 0$)

#### Demo 3: Propane at $320\text{ K}$, $1.0\text{ MPa}$ ($10\text{ bar}$) — Three Physical Roots & Phase Preference
- **State Classification:** Subcritical two-phase envelope ($\Delta < 0$, 3 physical roots $Z > B = 0.021165$)
- **Root 1 (Liquid-like):**
  - $Z = 0.036453$, $v = 9.6988 \times 10^{-5}\text{ m}^3/\text{mol}$ ($0.09699\text{ L/mol}$)
  - $\phi = 1.239280$, $f = 1.239280\text{ MPa}$
  - $H^{\text{res}} = -14\,854.34\text{ J/mol}$, $S^{\text{res}} = -48.2035\text{ J}/(\text{mol}\cdot\text{K})$
  - Mechanical Stability: STABLE ($(\partial P/\partial v)_T < 0$)
- **Root 2 (Intermediate):**
  - $Z = 0.091826$, $v = 2.4432 \times 10^{-4}\text{ m}^3/\text{mol}$
  - $\phi = 1.385615$, $H^{\text{res}} = -8210.60\text{ J/mol}$
  - Mechanical Stability: **UNSTABLE** ($(\partial P/\partial v)_T = +1.8825 \times 10^{10} > 0$); automatically discarded as unphysical.
- **Root 3 (Vapor-like):**
  - $Z = 0.850555$, $v = 2.2630 \times 10^{-3}\text{ m}^3/\text{mol}$ ($2.263\text{ L/mol}$)
  - $\phi = 0.868006$, $f = 0.868006\text{ MPa}$
  - $H^{\text{res}} = -1131.04\text{ J/mol}$, $S^{\text{res}} = -2.3575\text{ J}/(\text{mol}\cdot\text{K})$
  - Mechanical Stability: STABLE ($(\partial P/\partial v)_T < 0$)
- **Phase Equilibrium Analysis:**
  - $\Delta G^{\text{res}} = G_L - G_V = R T \ln(\phi_L / \phi_V) = +947.42\text{ J/mol} > 0$
  - Since $\phi_V < \phi_L$, the **Vapor phase has lower chemical potential and is the thermodynamically preferred stable state**.

---

## 6. Physical Assumptions & Known Model Limitations

1. **Pure Components Only:** Multicomponent gas/liquid mixtures and mixing rules (e.g. van der Waals one-fluid mixing rules) are not implemented.
2. **Canonical Peng–Robinson 1976 Formulation:** Uses the original 1976 acentric factor correlation $\kappa(\omega)$. Modifications for polar fluids or heavy hydrocarbons (e.g. Stryjek–Vera PRSV) are out of scope.
3. **Steam Reference State Limitation:** Gas-phase NIST Shomate parameters for Water Vapor ($\text{H}_2\text{O}$) are defined only on $[500.0, 1700.0]\text{ K}$. Because the standard reference state $T_{\text{ref}} = 298.15\text{ K}$ lies below $500\text{ K}$, total enthalpy and total entropy are reported as `N/A` with an explicit diagnostic message. PVT, fugacity, and residual departure properties remain fully accessible and exact.
4. **No Extrapolation:** Evaluating outside the NIST Shomate validity range rejects total property calculation via exceptions rather than silently extrapolating unphysical polynomials.
5. **Cubic Equation Density Limitations:** Like all two-parameter cubic equations of state, liquid molar volumes near the critical point exhibit typical deviations ($\sim 5\text{--}10\%$) without volume translation (Peneloux).

---

## 7. Project Structure
```text
thermo-eos-solver/
├── CMakeLists.txt              # Top-level CMake configuration
├── README.md                   # Project documentation & user guide
├── include/
│   ├── cli/
│   │   └── CLIController.hpp   # Interactive menu controller & formatters
│   └── thermo/
│       ├── Component.hpp       # PureComponent and Shomate structures
│       ├── ComponentDatabase.hpp # Database of 8 curated fluids
│       ├── Constants.hpp       # Physical constants & numerical tolerances
│       ├── CubicSolver.hpp     # Analytical Cardano/Vieta cubic solver
│       ├── PengRobinsonEOS.hpp # PR-1976 EOS engine & thermodynamic functions
│       └── StateProperties.hpp # PVTResult, PVTRoot, and phase classifications
├── src/
│   ├── cli_main.cpp            # Main entry point for thermo_cli binary
│   ├── cli/
│   │   └── CLIController.cpp   # Implementation of interactive CLI
│   ├── Component.cpp           # Component validation & Shomate integration
│   ├── ComponentDatabase.cpp   # Database initialization & lookup
│   ├── CubicSolver.cpp         # Cubic root solver implementation
│   └── PengRobinsonEOS.cpp     # PR-1976 thermodynamics implementation
├── tests/
│   ├── CMakeLists.txt          # CTest test target definitions
│   ├── test_cubic_solver.cpp   # Unit tests for cubic solver (48 assertions)
│   ├── test_components.cpp     # Unit tests for database & Shomate (44 assertions)
│   ├── test_pvt_eos.cpp        # Unit tests for PVT behavior (119 assertions)
│   ├── test_thermo_properties.cpp # Unit tests for H/S & fugacity (135 assertions)
│   └── test_cli.cpp            # Unit tests for CLI interface (27 assertions)
└── docs/                       # Comprehensive session audit reports
```
