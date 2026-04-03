# Microtopography-Aware Peatland Hydrology in NoahMP

## 1. Introduction

Land surface models typically represent the soil surface as a flat, horizontally uniform slab. In peatlands, however, the surface exhibits pronounced microtopography — alternating hummocks and hollows with elevation variations on the order of decimetres. This microtopography controls the partitioning of water between the soil matrix and free surface water ponded in hollows, the effective soil moisture seen by remote sensing, and the specific yield that governs water table dynamics.

PEAT-CLSM (Bechtold et al., 2019) introduced peatland-specific treatments of surface water storage and flooded fraction into the NASA Catchment Land Surface Model using empirical equations fitted to numerical simulations. Those relationships parameterised flooded fraction, soil-flux fraction ($f_\mathrm{soil}$), and surface water storage as functions of water table depth but lacked a closed-form physical derivation.

Here we describe a physically based peatland hydrology scheme in NoahMP. Following Dettmann & Bechtold (2015), all relationships are derived analytically from the integration of a standard soil water retention curve over a Gaussian microtopographic surface distribution. The scheme operates in two regimes — a hydrostatic-equilibrium path for shallow water tables and a Richards-equation path for deeper water tables — coupled through a head-shift mapping that links the one-dimensional flat-surface Richards domain to the microtopography-averaged column.

## 2. Theory

### 2.1 Microtopographic surface distribution

The surface elevation $z_s$ (positive upward from the mean surface) is modelled as a truncated Gaussian:

$$z_s \sim \mathcal{N}(0,\,\sigma^2), \quad z_s \in [-z_\mathrm{trunc},\, +z_\mathrm{trunc}]$$

with $\sigma = 0.16$ m and $z_\mathrm{trunc} = 1.0$ m. The cumulative distribution function (CDF)

$$F_s(z) = \Phi\!\left(\frac{z}{\sigma}\right)$$

gives the fraction of the surface area below elevation $z$, where $\Phi$ is the standard normal CDF.

### 2.2 Campbell soil water retention

Soil water retention follows the Campbell (1974) model:

$$\theta(h) = \begin{cases} \theta_s & h \ge -h_e \\ \theta_s \left(\dfrac{|h|}{h_e}\right)^{-1/b} & h < -h_e \end{cases}$$

where $h$ is the pressure head, $h_e$ is the air-entry suction, and $b$ is the pore-size distribution exponent. Under hydrostatic equilibrium with the water table at $z_\mathrm{wt}$, the pressure head at elevation $z$ is $h = z_\mathrm{wt} - z$. Campbell hydraulic conductivity follows $K(\theta) = K_s\,(\theta/\theta_s)^{2b+3}$.

### 2.3 Flooded fraction

The flooded fraction — the area fraction of the surface below the water table — follows directly from the surface CDF:

$$\alpha(z_\mathrm{wt}) = F_s(z_\mathrm{wt})$$

### 2.4 Surface water storage

The surface water storage per unit area [m] equals the volume of water ponded in hollows (Dettmann & Bechtold, 2015, Eq. 3):

$$W_\mathrm{surface}(z_\mathrm{wt}) = \int_{-z_\mathrm{trunc}}^{\min(z_\mathrm{wt},\, z_\mathrm{trunc})} F_s(z)\,dz$$

### 2.5 Soil water storage

Because the soil column extends to different depths across the microtopographic distribution, the effective soil water storage [m] integrates the retention curve weighted by the local soil fraction (Dettmann & Bechtold, 2015, Eq. 4):

$$W_\mathrm{soil}(z_\mathrm{wt}) = \int_{z_\mathrm{col,bot}}^{z_\mathrm{trunc}} \bigl[1 - F_s(z)\bigr]\;\theta(z_\mathrm{wt} - z)\;dz$$

At elevation $z$, the fraction $1 - F_s(z)$ is occupied by soil and the remainder by open water (when $z < z_\mathrm{wt}$).

### 2.6 Specific yield decomposition and flux partitioning

Dettmann & Bechtold (2015, Eq. 2) decompose the total specific yield into a soil and a surface component:

$$S_y = S_{y,\mathrm{soil}} + S_{y,\mathrm{surface}}$$

The soil component captures water released by desaturation of the soil matrix when the water table drops; the surface component captures the change in ponded water volume:

$$S_{y,\mathrm{surface}} = \frac{1}{\Delta z_\mathrm{wt}} \int_{z_l}^{z_u} F_s(z)\,dz$$

The fraction of vertical fluxes (infiltration, evaporation, transpiration) directed to the soil matrix is

$$f_\mathrm{soil} = \frac{S_{y,\mathrm{soil}}}{S_{y,\mathrm{soil}} + S_{y,\mathrm{surface}}}$$

evaluated over a small water table increment ($\pm 5$ mm) around the current depth. When the water table is well below the microtopographic range, $f_\mathrm{soil} \to 1$; as the water table rises and hollows flood, $f_\mathrm{soil}$ decreases because an increasing fraction of water table fluctuation is accommodated by changes in surface water storage rather than soil moisture.

All infiltration, evaporation, transpiration, and runoff fluxes are split accordingly: a fraction $f_\mathrm{soil}$ enters the soil Richards equation, and the complement $1 - f_\mathrm{soil}$ is routed to (or from) the surface water pool.

### 2.7 Numerical integration

All integrals over the microtopography distribution are evaluated with Gauss–Legendre quadrature on $[-z_\mathrm{trunc},\, z_\mathrm{trunc}]$. The standard normal CDF is computed via the Abramowitz & Stegun (1964, Eq. 7.1.26) rational approximation to the error function.

## 3. Two-Regime Algorithm

The soil water state is advanced in time with a two-regime algorithm that selects a hydrostatic-equilibrium path or a Richards-equation path depending on the water table depth at the beginning of each soil timestep.

### 3.1 Equilibrium path (WTD < 0.5 m)

When the water table is shallow, the soil column is close to hydrostatic equilibrium and capillary adjustment is fast relative to the model timestep. The equilibrium path bypasses the Richards equation and instead advances the soil state by direct inversion of the storage–water-table relationship.

The total water (soil + surface) is updated by the net flux balance:

$$W_\mathrm{total}^{n+1} = W_\mathrm{total}^n + (\text{Infiltration} - \text{Evaporation} - \text{Transpiration} - \text{Runoff}_\mathrm{sub})\,\Delta t$$

The new water table $z_\mathrm{wt}^{n+1}$ is found by numerically inverting $W_\mathrm{soil}(z_\mathrm{wt}) + W_\mathrm{surface}(z_\mathrm{wt}) = W_\mathrm{total}^{n+1}$ using a safeguarded Newton–Raphson method. Layer soil moisture is then set to the equilibrium profile at $z_\mathrm{wt}^{n+1}$, normalised to match the conserved soil water storage exactly. This path provides exact mass conservation without numerical diffusion.

### 3.2 Richards path (WTD $\ge$ 0.5 m)

When the water table is deeper, vertical redistribution may be far from equilibrium — infiltration wetting fronts, evaporative surface drying, and gravity drainage create profiles that deviate significantly from the hydrostatic shape. The Richards equation is solved explicitly in a flat-surface one-dimensional domain.

#### 3.2.1 Forward transfer (microtopography → flat)

Before the Richards solve, the column-averaged (microtopographic) soil moisture profile is translated into a flat-surface equivalent. For each active layer, a uniform pressure-head anomaly $\Delta h_i$ is diagnosed such that

$$\theta_\mathrm{micro,\,layer\,}i = \theta_\mathrm{eq,micro}(z_\mathrm{wt} - \Delta h_i)$$

i.e., the layer's microtopography-averaged moisture corresponds to shifting the equilibrium water table by $\Delta h_i$. The flat-domain point value at the layer midpoint $z_\mathrm{mid}$ is then computed from the shifted pressure head $h = z_\mathrm{wt} - \Delta h_i - z_\mathrm{mid}$ via the Campbell curve.

#### 3.2.2 Active/inactive layer classification

A layer is classified as *inactive* (decoupled from Richards) when the water table lies above the layer bottom. The topmost inactive layer index, $\mathrm{SatTopInd}$, is determined by

$$|z_\mathrm{bot,\,layer}| \ge \mathrm{WTD}$$

Inactive layers remain at saturation with zero inter-layer fluxes. The transitional layer (immediately above $\mathrm{SatTopInd}$) has its bottom flux removed to prevent spurious drainage across the water table.

#### 3.2.3 Richards equation with peatland modifications

The Richards equation is solved in head-gradient form:

$$\frac{\partial\theta}{\partial t} = \frac{\partial}{\partial z}\!\left[K(h)\left(\frac{\partial h}{\partial z} + 1\right)\right] - f_\mathrm{soil}\,S_\mathrm{sink}$$

where the sink term $S_\mathrm{sink}$ includes transpiration and soil evaporation. The head-gradient formulation (rather than the diffusivity form) provides better numerical stability in near-saturated peat soils. Bottom drainage is set to zero for peatlands. Only the $f_\mathrm{soil}$ fraction of transpiration and evaporation enters the Richards domain; the complement is accounted for in the surface water balance.

#### 3.2.4 Transpiration distribution

Transpiration is distributed across active soil layers following the original root-fraction weighting computed by NoahMP, but restricted to active layers only (layers above SatTopInd). The transpiration of inactive layers is redistributed to the active layers by scaling them proportionally, preserving the total transpiration amount.

#### 3.2.5 Subsurface runoff removal

After the Richards solve, the $f_\mathrm{soil}$ fraction of subsurface runoff is removed from active layers weighted by the product of layer hydraulic conductivity and thickness, $K_i\,\Delta z_i$, following the standard NoahMP approach but restricted to active layers.

#### 3.2.6 Backward transfer (flat → microtopography)

After the Richards solve and runoff removal, the flat-domain profile must be mapped back to the microtopographic column. The procedure is:

1. A conserved soil water target $W_\mathrm{target}$ is computed by applying the Richards-induced changes (deltas between pre- and post-Richards flat profiles) to the original microtopographic profile.
2. A diagnostic water table $\mathrm{WTD_{micro}}$ is found by inverting $W_\mathrm{soil}(z_\mathrm{wt}) = W_\mathrm{target}$.
3. Head-shift anomalies $\Delta h_i^\mathrm{flat}$ are extracted from the post-Richards profile at $\mathrm{WTD_{micro}}$.
4. A reference water table $\mathrm{WTD_{ref}}$ is found by bisection such that the microtopography-aware profile with the transferred head shifts reproduces $W_\mathrm{target}$ exactly.
5. The final layer moisture values are set to $\theta_\mathrm{micro}(z_\mathrm{wt,ref},\,\Delta h_i^\mathrm{flat})$, normalised for exact mass conservation.

This two-step mapping (forward + backward) allows the Richards equation to operate in its natural flat-surface domain while preserving microtopography-aware storage and water table dynamics.

### 3.3 Surface water balance

The change in surface water storage over a timestep is accumulated from the $(1 - f_\mathrm{soil})$ fractions of all fluxes:

$$\Delta W_\mathrm{surface} = (1 - f_\mathrm{soil})\,(\text{Infiltration} - \text{Evaporation} - \text{Transpiration} - \text{Runoff}_\mathrm{sub})\,\Delta t$$

The flooded fraction is updated from the final water table depth using the surface CDF. Total water — soil moisture plus surface water storage — is conserved by construction.

## References

- Abramowitz, M. & Stegun, I. A. (1964). *Handbook of Mathematical Functions*. Dover.
- Bechtold, M., De Lannoy, G. J. M., Koster, R. D., Reichle, R. H., Mahanama, S. P., Bleuten, W., Bourgault, M. A., Brümmer, C., Burdun, I., Desai, A. R., Devber, K., Gavazzi, M. J., Glagolev, M. V., Mezbahuddin, M., Ringgaard, R. & Humphreys, E. R. (2019). PEAT-CLSM: A specific treatment of peatland hydrology in the NASA Catchment Land Surface Model. *Journal of Advances in Modeling Earth Systems*, 11, 2130–2162.
- Campbell, G. S. (1974). A simple method for determining unsaturated conductivity from moisture retention data. *Soil Science*, 117, 311–314.
- Dettmann, U. & Bechtold, M. (2015). Deriving effective soil water retention characteristics from shallow water table fluctuations in peatlands. *Hydrological Processes*, 29, 3925–3940.
