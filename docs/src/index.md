
# AmphiDEB.jl

## Installation

AmphiDEB depends on `EcotoxSystems.jl`, which is curently not registered. 
Installation thus requires

```Julia
using Pkg
Pkg.add(url="htts://github.com/simonhansul/EcotoxSystems.jl")
Pkg.add(url="htts://github.com/simonhansul/AmphiDEB")
```

Guided examples are given in the [examples directory](https://github.com/SimonHansul/AmphiDEB/tree/main/examples).


## Notes on notation and terminology

### Simplified DEB notation

For DEB-TKTD parameters and state variables, 
we use a slightly simplified notation, following these rules:

- Leading capitalized character indicates a state variable
    - `X`: External resource
    - `X_emb`: Embryonic buffer
    - `I`: Ingested resource
    - `A`: Assimilates
    - `M`: Somatic maintenance costs
    - `H`: Maturity
    - `J`: Maturity maintenance costs
    - `R`: Reproduction buffer
- `k_X` is a rate constant related to state variable `X`
- `eta_XY` is the efficiency of transforming `X` into `Y`
- `X_int` is the initial value of `X`
- `y_X` indicates the response to `X` (e.g. chemical stressor or temperature). The relative response is *always* defined as a value that can be multiplied with the relevant term in the derivatives without further transformation.

Furthermore, surface areas and volumes are only given up to a proportionality:
- surface area $\propto$ mass $^{2/3}$
- volume $\propto$ mass

The zoom factor is based on the maximum structural masses of individuals.

### Life stage specificity

Some parameters are life stage specific, indicated by suffixes `_emb`, `_lrv`, `_mt`, `_juv`, `_ad`. <br>
The suffix indicates the life stage for which a switch in the value occurrs. <br>
For example, we have values `k_M_emb` and `k_M_juv`, which means that `k_M` has the value `k_M_emb` for embryos, larvae and metamorphs, 
and value `k_M_juv` for juveniles and adults.

## Input parameters

Listed below are the input parameters of the model (valid for v0.3.3). <br> 
For guidance on how to use these, please consult the examples. <br>
Please note that the default parameters are largely arbitrary and do not reflect a particular species.




### Global (`glb`)

| Parameter                  | Default value                    | Comment                                                                 |
|---------------------------|--------------------------|-------------------------------------------------------------------------|
| t_max                     | `56.`                    | Maximum simulation time simulation                                                       |
| N0                        | `1.`                     | Number of starting individuals
| dX_in                     | `[20., 20.]`             | Food input rate [mg d⁻¹], assuming *ad libitum* feeding                |
| k_V                       | `[0., 0.]`               | Dilution rate in the aquatic medium[d⁻¹]                                |
| V_patch                   | `[1., 1.]`               | Simulated volume and surface area of the aquatic and terrestic habitat, respectively
| T                         | `293.15`                 | Ambient temperature [K]                                                |
| C_W                       | `[0.;]`                  | Exposure concentrations (matrix: treatments × compounds)               |
| pathogen_inoculation_dose| `0.`                     | Pathogen spores added to medium [# spores]                              |
| pathogen_inoculation_time| `30.`                    | Time-point of pathogen inoculation [t]                                 |
| medium_renewals           | `[0.]`                   | Time-points for media renewals (removes spores)                         |


### Species-specific (`spc`)

| Parameter          | Default value                            | Comment                                                                                  |
|-------------------|----------------------------------|------------------------------------------------------------------------------------------|
| **Metaparameters** |                                  |                                                                                          |
| Z                 | `Dirac(1.)`                      | zoom factor, based on maximum structural masses                                          |
| propagate.zoom.dI_max_emb        | `1/3`                            | dI ∝ Z^(1/3)                                                                             |
| propagate.zoom.dI_max_lrv        | `1/3`                            |                                                                                          |
| propagate.zoom.dI_max_juv        | `1/3`                            |                                                                                          |
| propagate.zoom.X_emb_int         | `1.`                             | X_emb ∝ Z                                                                                |
| propagate.zoom.H_j1              | `1.`                             | H ∝ Z                                                                                    |
| propagate.zoom.H_p               | `1.`                             |                                                                                          |
| propagate.zoom.K_X_lrv           | `1/3`                            | K_X ∝ dI ∝ Z^(1/3)                                                                       |
| propagate.zoom.K_X_juv           | `1/3`                            |                                                                                          |
| **Physiological baseline (DEB)** |                   |                                                                                          |
| X_emb_int         | `1`                              | Initial vitellus (≈ dry mass of an egg)                                                  |
| K_X_lrv           | `1.`                             | Larval half-saturation constant for food uptake                                          |
| K_X_juv           | `1.`                             | Juvenile/adult half-saturation constant for food uptake                                  |
| dI_max_emb        | `1`                              | Embryonic maximum specific ingestion rate                                                |
| dI_max_lrv        | `1`                              | Larval maximum specific ingestion rate                                                   |
| dI_max_juv        | `1`                              | Juvenile/adult maximum specific ingestion rate                                           |
| kappa_emb         | `0.8`                            | Embryonic to metamorph allocation fraction                                               |
| kappa_juv         | `0.8`                            | Juvenile/adult allocation fraction to soma                                               |
| gamma             | `0.5`                            | Larval allocation to metamorphic reserves                                                |
| eta_IA            | `0.54`                           | Assimilation efficiency                                                                  |
| eta_AS_emb        | `0.4`                            | Embryonic growth efficiency                                                              |
| eta_AS_juv        | `0.4`                            | Juvenile/adult growth efficiency                                                         |
| eta_AR            | `0.95`                           | Reproduction efficiency                                                                  |
| eta_SA            | `0.8`                            | Shrinking efficiency                                                                     |
| k_M_emb           | `0.11`                           | Embryonic somatic maintenance rate                                                       |
| k_M_juv           | `0.11`                           | Juvenile/adult somatic maintenance rate                                                  |
| k_J_emb           | `0.027`                          | Embryonic maturity maintenance rate                                                      |
| k_J_juv           | `0.027`                          | Juvenile/adult maturity maintenance rate                                                 |
| H_j1              | `1`                              | Maturity at start of metamorphosis                                                       |
| H_p               | `55.`                            | Maturity at puberty                                                                      |
| T_A               | `8000.`                          | Arrhenius temperature (K)                                                                |
| T_ref             | `293.15`                         | Reference temperature                                                                    |
| b_T               | `40.`                            | Temperature effect on resource allocation (cf. [Romoli et al. (2025)](https://www.sciencedirect.com/science/article/pii/S0304380024003247))                                          |
| **TKTD Parameters** |                                 |                                                                                          |
| h_b               | `0.`                             | Background mortality                                                                     |
| KD                | `[0. 0. 0. 0. 0. 0. 0.;]`         | k_D per PMoA and stressor                                                                |
| B                 | `[2. 2. 2. 2. 2. 2. 2.;]`         | Slope parameters                                                                         |
| E                 | `[1e10 1e10 1e10 1e10 1e10 1e10 1e10;]` | Sensitivity thresholds                                                            |
| KD_h              | `[0.;]`                          | k_D for GUTS-SD module per stressor  |
| E_h               | `[1e10;]`                        | Sensitivity for GUTS-SD                                                                  |
| B_h               | `[1.;]`                          | Slope for GUTS-SD                                                                        |
| C_h               | `[1.;]`                          | Convert response to hazard rate                                                          |
| S_rel_crit        | `0.66`                           | Critical relative body mass loss                                                         |
| h_S               | `0.6`                            | Hazard rate below critical mass                                                          |
| a_max             | `truncated(Normal(15*365, 1.5*365), 0, Inf)` | Maximum age [d]                                                               |
| tau_R             | `365.`                           | Reproduction period [d]                                                                  |
| **Pathogen effects** |                              |                                                                                          |
| Chi               | `LogNormal(log(1)+1^2, 1)`       | Killing rate modifier (log-normal with mode 1, σ=1)                                      |
| E_P               | `[Inf, Inf, Inf, Inf]`           | Sensitivity thresholds for pathogen effects                                              |
| B_P               | `[2., 2., 2., 2.]`               | Slope parameters for pathogen effects                                                    |


### Pathogens (`pth`)

The model is adopted from [Drawert et al. (2017)](https://royalsocietypublishing.org/doi/full/10.1098/rsif.2017.0480). Default values are the averages of the suggested ranges.


| Parameter | Default value                              | Comment                                       |
|-----------|------------------------------------|-----------------------------------------------|
| gamma     | `geomean([1e-6, 1.])`              | Zoospore encounter rate [`spores individual⁻¹ d⁻¹`]                       |
| eta       | `harmmean([5., 20.])`              | Zoospore production rate                      |
| v0        | `0.5`                              | Zoospore encystment rate                      |
| f         | `0.5`                              | Host reinfection fraction                     |
| sigma0    | `harmmean([0.1, 0.5])`             | Sporangia killing rate                        |
| sigma1    | `harmmean([0.1, 0.5])`             | Density-dependent killing rate                |
| mu        | `harmmean([0.01, 1.5])`            | Zoospore death rate                           |


## API 

```@autodocs
Modules = [AmphiDEB]
Order   = [:type, :function]
```