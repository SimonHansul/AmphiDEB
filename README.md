# AmphiDEB.jl: Modelling amphibian and individual populations with Dynamic Energy Budgets


[![CI](https://github.com/SimonHansul/AmphiDEB.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/SimonHansul/AmphiDEB.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/SimonHansul/AmphiDEB.jl/graph/badge.svg?token=BL1CFR86M6)](https://codecov.io/gh/SimonHansul/AmphiDEB.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)


This package implements a Dynamic Energy Budget (DEB) model for amphibians based on DEBkiss. <br>
The implementation basis for this is [EcotoxSystems.jl](https://github.com/simonhansul/ecotoxsystems.jl).


## Changelog 


### v0.1.3

- Changed all life stage transitions to sigmoid functions and tweaked betas to optimize performance

### v0.1.4

- Updated dependencies
- Adjusted slopes for life stage transitions
- Some optimizations in TKTD submodel (can be improved further)
- Worked towards proper IBM tests


### v0.1.5-0.1.6

- Various small bugfixes and improved tests

### v0.1.7

- Fixed a typo in calculation of metamorphic reserve
- Added definition of alternative model version M2, which views metamorphic reserve as sub-compartment of structure (=>`E_mt` is subject to maintenance costs and contributes to surface area scaling)


### v0.1.8 

- Testing a slightly changed formulation for the default model

### v0.1.9

- Implemented PMoA 6: Decrease in kappa
- Adjusted default params to be compatible with log-logistic responses


### v0.2.0


**Breaking changes**:

- `glb.dX_in`, `glb.k_V` and `glb.V_patch` have to be provided as 2-element Vectors. The elements represent the aquatic and terrestric environment, respectively. 


**Non-breaking changes**

- Optimized model implementation.

### v0.2.1

- Added unit tests for temperature and pathogen effects
- Simplified pathogen effects by assuming response to pathogen to depend on `P_S` only, instead of `P_S/(S^(2/3))`. The latter kept causeing numerical instabilities with solvers other than Euler. TBD.

### v0.2.2

- Technical bugfixes, resolved merge conflices caused during previous update


### v 0.2.3 

- Default increasing log-logistic response is changed from $1/LL2(x,p)$ to $1-log(LL2(x, p))$. Results in more plausible y-values for common value of slope $b$, if the same priors are used for all PMoAs.
- Type of drc model used for sublethal effects can be changed through paramter `spc.drcmodel_sublethal` (1 = log-logistic, 2  = linear with threshold). 


### v 0.3.0

**Breaking changes**

- In the default derivatives, `calc_dE_mt_lrv` has been changed so that `eta_AS_lrv` and `y_G` also affects `dE_mt`. This leads to more patterns in the TKTD simulations when the PMoA is `G` (decrease in `eta_AS`).
- `spc.drcmodel_sublethal` is not a parameter anymore. Different TKTD model configuration are provided through modular functions
    - `AmphiDEB_ODE_with_loglogistic_TD!`
    - `AmphiDEB_ODE_with_linear_TD!`
    - `AmphiDEB_individual_ODE_with_loglogistic_TD!`
    - `AmphiDEB_individual_ODE_with_linear_TD!`
    - `TKTD_LL2!`
    - `TKTD_linear!`

**Non-breaking changes**

- Implemented additional PMoA. PMoAs are now (in this order):
    1. Decrease in growth efficiency
    2. Increase in maintenance costs (somatic + maturity)
    3. Decrease in assimilation efficiency
    4. Decrease in reproduction efficiency
    5. Decrease in maturity threshold for metamorphosis
    6. Increase in maturity threshold for metamorphosis
    7. Decrease in $\kappa$

### v0.3.1

- Connection to `C_W` is now life stage specific: Only larvae are affected by aquatic exposure.

### v0.3.2

- bugfix in `TK_minimal_aquatic`

### v0.3.3

- LL2pos has a caveat to catch negative x values