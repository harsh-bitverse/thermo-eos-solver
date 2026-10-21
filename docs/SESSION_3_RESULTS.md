# CHL2002 Thermodynamics Project — Session 3: Peng–Robinson EOS Engine and PVT Calculations

**Document Version:** 1.0  
**Date:** October 21, 2026  
**Status:** Completed & Validated Module Report  

---

## 1. Executive Summary

In Session 3, the core thermodynamic engine for the pure-component Peng–Robinson (PR-1976) Equation of State was implemented in C++17. The engine provides:
1. **Analytical EOS Parameters:** Exact precomputation of critical parameters ($a_c, b$), acentric parameter ($\kappa$), reduced temperature ($T_r$), temperature-dependent alpha function ($\alpha(T_r)$), attractive parameter ($a(T)$), and analytical temperature derivative ($da/dT$).
2. **Dimensionless Polynomial Construction:** Accurate translation to $A(T, P)$ and $B(T, P)$, constructing the monic PR cubic polynomial in $Z$.
3. **Forward PVT Engine:** Solves the cubic equation using the existing `CubicSolver`, filters for physically admissible roots ($Z > B \iff v > b$), calculates molar volume $v = ZRT/P$, analytically computes the isothermal derivative $(dP/dv)_T$, and classifies roots by mechanical stability and phase behavior.
4. **Independent Inverse PVT Engine:** Directly evaluates pressure $P(T, v) = \frac{RT}{v-b} - \frac{a(T)}{v^2 + 2bv - b^2}$ with domain validation ($v > b$).
5. **Comprehensive Verification:** 113 automated unit test assertions in `tests/test_pvt_eos.cpp` validating parameter dimensions, forward/inverse reciprocity ($|P_{\text{calc}} - P|/P \le 10^{-10}$), ideal-gas asymptotic limits ($P \to 0$), audited literature benchmarks (Methane, $\text{CO}_2$, Propane), derivative cross-validation, and domain guardrails.

---

## 2. Governing Equations & Formulations

### 2.1 Pure-Component Peng–Robinson Parameters (1976)
$$P = \frac{RT}{v - b} - \frac{a(T)}{v(v + b) + b(v - b)} = \frac{RT}{v - b} - \frac{a(T)}{v^2 + 2bv - b^2}$$

* **Universal Gas Constant:** $R = 8.31446261815324 \text{ J}\cdot\text{mol}^{-1}\cdot\text{K}^{-1}$ (CODATA 2018/2022).
* **Critical Constants:**
  $$a_c = \Omega_a \frac{R^2 T_c^2}{P_c}, \quad \Omega_a = 0.4572355289213822... \approx 0.45724$$
  $$b = \Omega_b \frac{R T_c}{P_c}, \quad \Omega_b = 0.07779607390388846... \approx 0.07780$$
  $$\kappa = 0.37464 + 1.54226\,\omega - 0.26992\,\omega^2$$
  $$T_r = \frac{T}{T_c}$$
  $$\alpha(T_r, \omega) = \left[ 1 + \kappa \left( 1 - \sqrt{T_r} \right) \right]^2$$
  $$a(T) = a_c \cdot \alpha(T_r, \omega)$$

### 2.2 Analytical Temperature Derivative $\frac{da}{dT}$
$$\frac{da}{dT} = - \frac{a_c \kappa \left[ 1 + \kappa \left( 1 - \sqrt{T_r} \right) \right]}{\sqrt{T \cdot T_c}} = - \frac{a_c \kappa \sqrt{\alpha(T_r)}}{\sqrt{T \cdot T_c}} \quad [\text{Pa}\cdot\text{m}^6\cdot\text{mol}^{-2}\cdot\text{K}^{-1}]$$

### 2.3 Dimensionless PR Cubic Polynomial
Substituting $v = \frac{Z R T}{P}$ with $A = \frac{a(T) P}{(R T)^2}$ and $B = \frac{b P}{R T}$:
$$Z^3 - (1 - B)Z^2 + (A - 3B^2 - 2B)Z - (AB - B^2 - B^3) = 0$$

* $c_2 = -(1 - B)$
* $c_1 = A - 3B^2 - 2B$
* $c_0 = -(AB - B^2 - B^3)$

### 2.4 Isothermal Derivative $\left(\frac{\partial P}{\partial v}\right)_T$ (Mechanical Stability)
$$\left(\frac{\partial P}{\partial v}\right)_T = -\frac{RT}{(v - b)^2} + \frac{a(T)(2v + 2b)}{(v^2 + 2bv - b^2)^2} \quad [\text{Pa}\cdot\text{mol}\cdot\text{m}^{-3}]$$
* **Stable Root:** $\left(\frac{\partial P}{\partial v}\right)_T < 0$ (compression causes pressure to rise).
* **Unstable Root:** $\left(\frac{\partial P}{\partial v}\right)_T \ge 0$ (spinodal unobservable transition state).

---

## 3. C++ Architecture and Public Interfaces

### 3.1 `StateProperties.hpp`
Defines the thermodynamic state containers:
```cpp
namespace thermo {

enum class PhaseTag { SINGLE_PHASE, VAPOR, LIQUID, INTERMEDIATE };
enum class MechanicalStability { STABLE, UNSTABLE };

struct PVTRoot {
    double Z{1.0};
    double molar_volume{0.0};
    double poly_residual{0.0};
    double dP_dv{0.0};
    MechanicalStability stability{MechanicalStability::STABLE};
    PhaseTag phase{PhaseTag::SINGLE_PHASE};
    bool is_physical{false};
};

struct PVTResult {
    double T{0.0};
    double P{0.0};
    double A{0.0};
    double B{0.0};
    double discriminant{0.0};
    std::vector<PVTRoot> roots;

    [[nodiscard]] bool has_multiple_roots() const noexcept;
    [[nodiscard]] const PVTRoot* vapor_root() const noexcept;
    [[nodiscard]] const PVTRoot* liquid_root() const noexcept;
    [[nodiscard]] const PVTRoot* intermediate_root() const noexcept;
};

struct EOSParameters {
    double a, b, da_dT, ac, kappa, alpha, Tr;
};

} // namespace thermo
```

### 3.2 `PengRobinsonEOS.hpp`
Defines the core EOS engine:
```cpp
namespace thermo {

class PengRobinsonEOS {
public:
    explicit PengRobinsonEOS(PureComponent component);

    [[nodiscard]] const PureComponent& component() const noexcept;
    [[nodiscard]] double ac() const noexcept;
    [[nodiscard]] double b() const noexcept;
    [[nodiscard]] double kappa() const noexcept;
    [[nodiscard]] double Tr(double T) const;
    [[nodiscard]] double alpha(double T) const;
    [[nodiscard]] double a(double T) const;
    [[nodiscard]] double da_dT(double T) const;
    [[nodiscard]] double A_parameter(double T, double P) const;
    [[nodiscard]] double B_parameter(double T, double P) const;
    [[nodiscard]] EOSParameters parameters(double T) const;

    // Forward PVT: (T, P) -> Z, v, stability, classification
    [[nodiscard]] PVTResult calculate_pvt(double T, double P) const;

    // Inverse PVT: (T, v) -> P
    [[nodiscard]] double calculate_pressure(double T, double v) const;

    // Isothermal derivative (dP/dv)_T
    [[nodiscard]] double calculate_dP_dv(double T, double v) const;
};

} // namespace thermo
```

---

## 4. Root Filtering & Classification Logic

1. **Mathematical Solving:** Monic cubic polynomial is solved via `CubicSolver::solve(c2, c1, c0)` using Vieta's trigonometric method ($\Delta < 0$) or Cardano's method ($\Delta \ge 0$).
2. **Physical Admissibility Filter:**
   * Molar volume must satisfy $v > b \iff Z > B$.
   * A small numerical safety margin ($1.0 \times 10^{-12}$) prevents boundary roundoff:
     $$Z > B + 1.0 \times 10^{-12}$$
   * Roots violating this condition are non-physical and discarded.
3. **Molar Volume & Mechanical Stability:**
   * For each physical root: $v = \frac{ZRT}{P}$.
   * Analytical derivative $(dP/dv)_T$ is computed.
   * If $(dP/dv)_T < 0 \implies \text{STABLE}$; else $\text{UNSTABLE}$.
4. **Phase Tagging:**
   * If 1 physical root: Tagged as `SINGLE_PHASE` (supercritical or outside spinodal).
   * If 2 physical roots: Smallest root tagged as `LIQUID`, largest as `VAPOR`.
   * If 3 physical roots ($Z_0 < Z_1 < Z_2$):
     * $Z_0$ (smallest): `LIQUID` (dense phase, $(dP/dv)_T < 0$).
     * $Z_1$ (middle): `INTERMEDIATE` (mechanically unstable, $(dP/dv)_T > 0$).
     * $Z_2$ (largest): `VAPOR` (expanded phase, $(dP/dv)_T < 0$).
   * *Important Design Rule:* Root classification is strictly distinguished from global thermodynamic phase stability. Global stability requires equal-fugacity or Gibbs energy minimization (Session 4).

---

## 5. Independent Benchmark Calculations

Calculated using pure-component PR-1976 equations ($R = 8.31446261815324\text{ J}/(\text{mol}\cdot\text{K})$):

| Case ID | Substance | State ($T, P$) | Dimensionless Parameters | Depressed Cubic & Discriminant | Roots & Properties ($Z, v, \partial P/\partial v$) | Classification & Stability |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **BM-01** | **Methane** | $300\text{ K}$<br>$5.0\text{ MPa}$ | $A = 0.16248591$<br>$B = 0.05372479$ | $p = -0.25210166$<br>$q = -0.05382471$<br>$\Delta = +1.3085 \times 10^{-4} > 0$ | $Z = 0.901844$<br>$v = 4.4990 \times 10^{-4}\text{ m}^3/\text{mol}$<br>$\frac{\partial P}{\partial v} = -1.015 \times 10^{10}\text{ Pa}/(\text{m}^3/\text{mol})$ | **SINGLE_PHASE**<br>Mechanically STABLE |
| **BM-02** | **$\text{CO}_2$** | $300\text{ K}$<br>$4.0\text{ MPa}$ | $A = 0.25725664$<br>$B = 0.04276376$ | $p = -0.13919084$<br>$q = -0.02102135$<br>$\Delta = +1.0597 \times 10^{-5} > 0$ | $Z = 0.752334$<br>$v = 4.6914 \times 10^{-4}\text{ m}^3/\text{mol}$<br>$\frac{\partial P}{\partial v} = -6.098 \times 10^9\text{ Pa}/(\text{m}^3/\text{mol})$ | **SINGLE_PHASE**<br>Mechanically STABLE |
| **BM-03** | **$\text{CO}_2$** | $300\text{ K}$<br>$7.0\text{ MPa}$ | $A = 0.45019912$<br>$B = 0.07483658$ | $p = -0.00158471$<br>$q = +0.00116819$<br>$\Delta = +3.4102 \times 10^{-7} > 0$ | $Z = 0.198057$<br>$v = 7.0575 \times 10^{-5}\text{ m}^3/\text{mol}$<br>$\frac{\partial P}{\partial v} = -8.802 \times 10^{10}\text{ Pa}/(\text{m}^3/\text{mol})$ | **SINGLE_PHASE**<br>(Compressed Liquid)<br>Mechanically STABLE |
| **BM-04** | **Propane** | $320\text{ K}$<br>$1.0\text{ MPa}$ | $A = 0.15613051$<br>$B = 0.02116534$ | $p = -0.20691652$<br>$q = -0.03562470$<br>$\Delta = -1.0832 \times 10^{-5} < 0$ | **Root 0:** $Z = 0.036453, v = 9.699 \times 10^{-5}, \partial P/\partial v < 0$<br>**Root 1:** $Z = 0.091826, v = 2.443 \times 10^{-4}, \partial P/\partial v > 0$<br>**Root 2:** $Z = 0.850555, v = 2.263 \times 10^{-3}, \partial P/\partial v < 0$ | **Root 0: LIQUID, STABLE**<br>**Root 1: INTERMEDIATE, UNSTABLE**<br>**Root 2: VAPOR, STABLE** |

---

## 6. Test Results & Numerical Tolerances

Automated testing was conducted via CTest and `test_pvt_eos.exe`:

* **Total Assertions in `test_pvt_eos`:** 113 passed, 0 failed.
* **Test A (Parameter Consistency):** Verified parameter dimensions, finite bounds, and analytical $da/dT < 0$ across all gases.
* **Test B (Reciprocity):** Across 12 distinct $(T, P)$ states spanning Methane, $\text{CO}_2$, and Propane (from $1\text{ bar}$ to $200\text{ bar}$), forward/inverse reciprocity satisfied:
  $$\frac{|P_{\text{calc}} - P|}{P} \le 1.0 \times 10^{-10}$$
* **Test C (Ideal-Gas Limit):** Evaluated from $P = 100\text{ Pa}$ down to $P = 10^{-4}\text{ Pa}$. Confirmed $Z \to 1.0$ and $v \to RT/P$ with smooth asymptotic convergence.
* **Test D (Topology & Benchmarks):** Verified exact matching against audited values (Case 1: $Z = 0.901844$; Case 2: $Z = 0.752334$; Case 3: $Z = 0.198057$; Case 4: 3 roots with liquid, intermediate, vapor tagging).
* **Test E (Mechanical Stability):** Verified analytical $(dP/dv)_T$ against central finite difference $\frac{P(v+\delta) - P(v-\delta)}{2\delta}$ ($\delta = 10^{-5} v$) within relative error $\le 1.0 \times 10^{-5}$. Verified positive sign for intermediate root and negative sign for liquid and vapor roots.
* **Test F (Invalid Inputs):** Verified that $T \le 0$, $P \le 0$, $v \le b$, and non-finite numbers (NaN, Inf) reliably throw `std::invalid_argument`.

---

## 7. Known Limitations

1. **Global Phase Stability:** Session 3 classifies roots based on density ($Z_{\min}$ vs $Z_{\max}$) and mechanical stability ($dP/dv < 0$), but intentionally does not declare which phase is globally thermodynamically stable in the two-phase region. Molar Gibbs free energy / fugacity coefficient comparison ($f_{\text{liq}}$ vs $f_{\text{vap}}$) is reserved for Session 4.
2. **Cubic Equation Domain:** Cubic EOS models describe classical fluid behavior; extremely high pressures ($P > 500\text{ MPa}$) where solid phases form are beyond the physical scope of the PR EOS.

---

## 8. Git Commit Record

* **Commit Hash:** Pending final Session 3 commit.
* **Files Tracked:**
  * `include/thermo/StateProperties.hpp`
  * `include/thermo/PengRobinsonEOS.hpp`
  * `src/PengRobinsonEOS.cpp`
  * `tests/test_pvt_eos.cpp`
  * `CMakeLists.txt`
  * `tests/CMakeLists.txt`
  * `docs/IMPLEMENTATION_ROADMAP.md`
  * `docs/SESSION_3_RESULTS.md`
