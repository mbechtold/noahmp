# Microtopography-Aware Soil Water Retention and Surface Water Storage for Peatlands in NoahMP

## 1. Introduction and Motivation

Land surface models typically represent the soil column as a flat, horizontally uniform slab. In peatlands, however, the surface exhibits pronounced microtopography — alternating hummocks and hollows with elevation variations on the order of decimeters. This microtopography exerts a first-order control on the partitioning of water between the soil matrix and free surface water ponded in hollows, on the effective soil moisture observed by remote sensing, and on the specific yield that governs water table fluctuations.

The PEATCLSM model (Bechtold et al., 2019) introduced peatland-specific treatments of surface water storage and flooded fraction into the NASA Catchment Land Surface Model using empirical equations fitted to numerical simulations. These relationships parameterised the flooded fraction, the soil-flux fraction (f\_soil), and the surface water storage as functions of water table depth, but lacked a closed-form physical derivation, making them difficult to transfer across sites with different microtopographic characteristics or soil properties.

Here we implement a physically based framework in NoahMP following Dettmann & Bechtold (2015, *Hydrological Processes*), who derived analytical expressions for specific yield, surface water storage, and effective soil water content over a Gaussian microtopographic surface under hydrostatic equilibrium. Compared to the empirical PEATCLSM approach, the present implementation offers three key advantages: (1) all relationships follow directly from the integration of a standard soil water retention curve over a measurable microtopographic distribution, providing full mathematical transparency; (2) the approach depends on only three measurable quantities — the standard deviation of surface elevations (σ), the Campbell soil hydraulic parameters (θ\_s, h\_e, b), and the water table depth — making it straightforward to apply at any peatland site; and (3) the specific yield decomposition into soil and surface components provides a physically consistent partitioning of vertical fluxes between the soil matrix and ponded surface water.

## 2. Theory

### 2.1 Microtopographic surface distribution

We model the surface elevation z\_s as a Gaussian random variable with zero mean (referenced to the mean surface) and standard deviation σ = 0.16 m, truncated at ±z\_trunc = ±1.0 m:

$$z_s \sim \mathcal{N}(0, \sigma^2), \quad z_s \in [-z_{\mathrm{trunc}}, +z_{\mathrm{trunc}}]$$

The cumulative distribution function (CDF) of surface elevations, F\_s(z), gives the fraction of the surface below elevation z:

$$F_s(z) = \Phi\!\left(\frac{z}{\sigma}\right), \quad z \in [-z_{\mathrm{trunc}}, z_{\mathrm{trunc}}]$$

where Φ is the standard normal CDF. F\_s is clamped to 0 for z ≤ −z\_trunc and to 1 for z ≥ z\_trunc.

### 2.2 Flooded fraction

The flooded fraction — the fraction of the peatland surface covered by free-standing water — equals the fraction of the surface whose elevation lies below the water table:

$$\alpha(z_{\mathrm{wt}}) = F_s(z_{\mathrm{wt}})$$

When the water table is deep (z\_wt ≪ 0), α → 0; when the water table is well above the mean surface (z\_wt ≫ 0), α → 1. In contrast to the empirical flooded-fraction relationship used in PEATCLSM, this expression follows directly from the assumed Gaussian surface distribution.

### 2.3 Soil hydraulic model

The soil water retention follows the Campbell (1974) model as used throughout NoahMP:

$$\theta(h) = \begin{cases} \theta_s & \text{if } h \geq -h_e \\ \theta_s \left(\frac{|h|}{h_e}\right)^{-1/b} & \text{if } h < -h_e \end{cases}$$

where h is the pressure head, h\_e is the (positive) air-entry suction head, and b is the Campbell pore-size distribution exponent. Under hydrostatic equilibrium with water table at z\_wt, the pressure head at elevation z is h = z\_wt − z.

For peat soil (NoahMP soil type 17), the parameters are θ\_s = 0.880, h\_e = 0.024 m, b = 7.4.

### 2.4 Surface water storage

Following Eq. 3 of Dettmann & Bechtold (2015), the surface water storage per unit horizontal area [m] equals the depth of water ponded in hollows:

$$V_{\mathrm{surface}}(z_{\mathrm{wt}}) = \int_{-z_{\mathrm{trunc}}}^{\min(z_{\mathrm{wt}}, z_{\mathrm{trunc}})} F_s(z) \, dz$$

For z\_wt above z\_trunc, an additional open-water depth (z\_wt − z\_trunc) is added for the fully flooded portion.

### 2.5 Soil water storage with microtopography

The soil water storage [m] is obtained by integrating the soil moisture over the heterogeneous soil column, weighting by the fraction of ground present at each elevation (Eq. 4 in D&B 2015):

$$A_{\mathrm{soil}}(z_{\mathrm{wt}}) = \int_{-z_{\mathrm{trunc}}}^{+z_{\mathrm{trunc}}} \bigl[1 - F_s(z)\bigr] \, \theta(z_{\mathrm{wt}} - z) \, dz$$

At elevation z, the fraction (1 − F\_s(z)) is occupied by soil and the fraction F\_s(z) may be open water (if z < z\_wt). The pressure head at elevation z is h = z\_wt − z under hydrostatic equilibrium.

### 2.6 Specific yield decomposition

Dettmann & Bechtold (2015, Eq. 2) decompose the total specific yield S\_y into a soil component and a surface (open-water) component:

$$S_y(z_{\mathrm{wt}}) = S_{y,\mathrm{soil}}(z_{\mathrm{wt}}) + S_{y,\mathrm{surface}}(z_{\mathrm{wt}})$$

For a finite water table change from z\_l to z\_u:

$$S_{y,\mathrm{soil}} = \frac{1}{z_u - z_l} \int_{-z_{\mathrm{trunc}}}^{+z_{\mathrm{trunc}}} [1 - F_s(z)] \, [\theta(z_u - z) - \theta(z_l - z)] \, dz$$

$$S_{y,\mathrm{surface}} = \frac{1}{z_u - z_l} \int_{z_l}^{z_u} F_s(z) \, dz$$

### 2.7 Flux partitioning (f\_soil)

The fraction of vertical fluxes (infiltration, evaporation, transpiration) directed to the soil matrix versus the surface water is determined by the local specific yield decomposition:

$$f_{\mathrm{soil}} = \frac{S_{y,\mathrm{soil}}}{S_{y,\mathrm{soil}} + S_{y,\mathrm{surface}}}$$

evaluated over a small increment (±5 mm) around the current water table. When the water table is deep (z\_wt ≪ −σ), f\_soil → 1 (all flux to soil). As the water table rises and hollows begin to flood, f\_soil decreases progressively because an increasing fraction of water table fluctuation is accommodated by changes in surface water storage rather than soil moisture. This transition begins well below the mean surface — as soon as the water table enters the microtopographic range — and f\_soil continues to decrease toward zero as the water table rises above the mean surface and most of the area is inundated.

### 2.8 Transient soil moisture and equilibrium diagnostics

A key design choice in this implementation is that the transient soil moisture profile produced by the Richards equation is preserved as the model state. The equilibrium assumption (hydrostatic profile above the water table) is used only internally to diagnose the water table depth from the total soil water deficit, and to compute f\_soil and the flooded fraction.

This ensures that the modelled surface soil moisture retains realistic transient features — infiltration wetting fronts, evaporative drying at the surface, gravity redistribution — that are important for comparison with satellite-derived soil moisture products. At deeper water tables, the equilibrium assumption for the unsaturated zone becomes less accurate, making the preservation of the Richards-solved profile particularly valuable.

For offline diagnostic purposes, the effective (microtopography-integrated) water content for a model layer can be computed as:

$$\theta_{\mathrm{eff}}^{[z_{\mathrm{bot}}, z_{\mathrm{top}}]} = \frac{1}{\Delta z} \int_{z_{\mathrm{bot}}}^{z_{\mathrm{top}}} \bigl[(1-F_s(z))\,\theta(z_{\mathrm{wt}}-z) + F_s(z)\,\mathbb{1}_{z < z_{\mathrm{wt}}}\bigr] \, dz$$

where the first term is the soil-matrix contribution and the second is the surface-water contribution in hollows below the water table. This function is available in `PeatMicroTopoMod` but is not applied to the prognostic state variables.

### 2.9 Numerical integration

All integrals over the microtopography distribution are evaluated using 20-point Gauss–Legendre quadrature on the truncated interval [−z\_trunc, +z\_trunc]. The standard normal CDF is computed via the Abramowitz & Stegun (1964, Eq. 7.1.26) rational approximation to the error function, with maximum error 1.5 × 10⁻⁷.

## 3. Implementation in NoahMP

### 3.1 Sign convention

NoahMP's `WaterTableDepth` is positive downward. Internally, we convert to the Dettmann & Bechtold (2015) convention where z is positive upward from the mean surface:

$$z_{\mathrm{wt}} = -\text{WaterTableDepth}$$

The allowed range is WaterTableDepth ∈ [−1.0, +3.0] m, corresponding to z\_wt ∈ [−3.0, +1.0] m (water table from 3 m below surface to 1 m above the highest hummock).

### 3.2 New module: PeatMicroTopoMod.F90

This self-contained module implements all core physics functions:

| Function | Description |
|----------|-------------|
| `erf_approx(x)` | Error function approximation (Abramowitz & Stegun 7.1.26) |
| `phi_normal(x)` | Standard normal CDF |
| `Fs_cdf(z)` | CDF of surface elevations |
| `FloodedFrac(z_wt)` | Flooded fraction α(z\_wt) |
| `theta_campbell(h, θ_s, h_e, b)` | Campbell retention curve |
| `SurfaceWaterStorage(z_wt)` | Surface water volume [m] |
| `SoilWaterStorageMicroTopo(z_wt, ...)` | Soil water storage [m] |
| `TotalWaterStorageMicroTopo(z_wt, ...)` | Total (soil + surface) storage [m] |
| `EffectiveSoilMoistureLayer(z_bot, z_top, z_wt, ...)` | Effective θ for a model layer |
| `SysoilMicroTopo(z_l, z_u, ...)` | Soil component of specific yield |
| `SysurfaceMicroTopo(z_l, z_u)` | Surface component of specific yield |
| `SytotalMicroTopo(z_l, z_u, ...)` | Total specific yield |
| `FsoilMicroTopo(z_wt, ...)` | Soil flux fraction f\_soil |
| `FindWaterTable(W_total, ...)` | Bisection for z\_wt from total water |
| `SoilDeficitMicroTopo(z_wt, ...)` | Soil water deficit for equilibrium |

Module parameters: σ = 0.16 m, z\_trunc = 1.0 m, 20-point Gauss–Legendre nodes and weights.

### 3.3 Modified module: MicroTopoCorrectionMod.F90

Computes `FloodedFraction` and `f_soil` using the physics-based functions from `PeatMicroTopoMod`. Called before the Richards equation solver at each soil timestep.

### 3.4 Modified module: WaterTableEquilibriumPeatMod.F90

Finds the equilibrium water table using a bisection method on the microtopography-aware soil deficit function `SoilDeficitMicroTopo`. With microtopography extending to +1 m above the mean surface, there is always unsaturated soil in the hummocks, so the deficit is strictly positive under any realistic condition and a single deficit-based bisection suffices.

### 3.5 Modified module: SoilWaterMainMod.F90

The main soil water orchestrator was updated with the following changes:

1. **Minimum water table depth** changed from −0.2449 m to −1.0 m (full truncation range)
2. **Pre-Richards (predictor):** `MicroTopoCorrection` computes f\_soil and FloodedFraction at the current water table depth; infiltration is split into f\_soil (to soil) and (1−f\_soil) (to FSW\_change); the water table depth at this point is saved as WTD\_begin
3. **Richards equation:** transpiration and evaporation are scaled by f\_soil (only the soil fraction of fluxes is processed)
4. **Post-Richards FSW\_change:** accumulates (1−f\_soil) fractions of infiltration, evaporation, transpiration, and runoff removal
5. **WTD diagnosis and corrector step:** After Richards and runoff removal, `WaterTableEquilibriumPeat` diagnoses a new WTD\_end from the updated soil moisture deficit. A corrector then recomputes f\_soil at the midpoint water table (WTD\_begin + WTD\_end)/2. The difference between the corrected and original f\_soil, multiplied by the total net flux, gives a flux correction [mm] that is redistributed between the soil column and surface water storage. WTD is then re-diagnosed from the corrected soil moisture. This predictor–corrector approach aligns the water level changes implied by the soil-side and surface-water-side flux partitions, reducing errors from the frozen f\_soil assumption without the cost of re-running Richards.
6. **Transient profile preserved:** The Richards-solved soil moisture profile is retained as the model state — no equilibrium overwrite is applied. The equilibrium assumption is used only internally for WTD diagnosis, f\_soil, and FloodedFraction.

### 3.6 Modified module: SoilWaterDiffusionRichardsMod.F90

For peatlands (OptPeatlandPhysics = 1), the transpiration and soil evaporation terms in the WaterExcess computation are multiplied by f\_soil, so only the soil-directed fraction of these fluxes enters the Richards equation. Bottom drainage is set to zero.

### 3.7 Modified module: RunoffSubSurfacePeatlandMod.F90

The Ivanov-based subsurface runoff scheme was updated to:
- Always call `WaterTableEquilibriumPeat` (removed the previous WTD > 0.1 threshold guard)
- Use the extended WTD range (clamp at −1.0 m instead of −0.2449 m) for transmissivity computation

### 3.8 New state variable: FloodedFraction

Added to `WaterVarType` as `FloodedFraction` (fraction of the surface below the water table), computed from the Gaussian CDF. The former `AR1` variable has been removed.

### 3.9 Water balance

The water balance tracker (`BalanceErrorCheckMod`) accounts for surface water storage changes through `FSW_change`, which accumulates the (1−f\_soil) fraction of all fluxes: infiltration, evaporation, transpiration, and runoff removal. The corrector step (Section 3.5) redistributes water between `FSW_change` and `SoilLiqWater` but conserves total water exactly: the amount subtracted from `FSW_change` equals the amount added to the soil column. The transient Richards profile in `SoilLiqWater` tracks the soil-side water, so `WaterStorageTotEnd` (which sums `SoilMoisture × thickness`) plus `FSW_change` closes the water balance.

## 4. Summary of Code Changes

| File | Status | Description |
|------|--------|-------------|
| PeatMicroTopoMod.F90 | **New** | Core physics: Gaussian CDF, Campbell retention, GL quadrature, all D&B 2015 functions |
| MicroTopoCorrectionMod.F90 | **Revised** | f\_soil and FloodedFraction from Sy decomposition |
| WaterTableEquilibriumPeatMod.F90 | **Revised** | Equilibrium WTD with microtopography deficit |
| SoilWaterMainMod.F90 | **Revised** | Flux partitioning with predictor–corrector, extended WTD range, transient SM preserved |
| SoilWaterDiffusionRichardsMod.F90 | **Revised** | f\_soil scaling of evaporation/transpiration |
| RunoffSubSurfacePeatlandMod.F90 | **Revised** | Extended WTD range for transmissivity |
| WaterVarType.F90 | **Revised** | Added FloodedFraction |
| WaterVarInitMod.F90 | **Revised** | Initialize FloodedFraction |
| Makefile | **Revised** | Added new .o files and dependencies |

## References

- Abramowitz, M., & Stegun, I. A. (1964). *Handbook of Mathematical Functions*. Dover.
- Bechtold, M., et al. (2019). PEAT-CLSM: A specific treatment of peatland hydrology in the NASA Catchment Land Surface Model. *Journal of Advances in Modeling Earth Systems*, 11, 2130–2162.
- Campbell, G. S. (1974). A simple method for determining unsaturated conductivity from moisture retention data. *Soil Science*, 117, 311–314.
- Dettmann, U., & Bechtold, M. (2015). Deriving effective soil water retention characteristics from shallow water table fluctuations in peatlands. *Hydrological Processes*, 29, 3925–3940. DOI: 10.1002/hyp.10475
