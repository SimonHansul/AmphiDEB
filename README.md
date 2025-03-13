# AmphiDEB.jl: Modelling amphibian and individual populations with Dynamic Energy Budgets


[![CI](https://github.com/SimonHansul/AmphiDEB.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/SimonHansul/AmphiDEB.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/SimonHansul/AmphiDEB.jl/graph/badge.svg?token=BL1CFR86M6)](https://codecov.io/gh/SimonHansul/AmphiDEB.jl)

## TKTD interface

The TKTD module is designed so that mixtures with an arbitrary number of chemicals can be simulated, 
and each chemical can act via arbitrary combinations of PMoAs. <br>
For this pupose, TKTD parameters for sublethal effects are stored in matrices, where the rows represent chemicals and the columns represent PMoAs. <br>
This affects the parameters $k_{D_j}$, $e_j$ and $b_j$. <br>
The PMoAs have a fixed order in the parameter matrices:

1. Decrease in growth efficiency ($G$)
2. Increase in somatic and maturity maintenance costs ($M$)
3. Decrease in assimilation efficiency ($A$)
4. Decrease in reproduction efficiency ($R$)
5. Increase in maturity threshold for metamorphosis ($H$)
6. Decrease in $\kappa$ ($\kappa$, acceleration of ontogenesis)

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

- Testing a slightly changed formulation for the default model (1-$\gamma$ pulled out of the parenthesis so that $\gamma$ does not effect Equilbirium $S$)

### v0.1.9

- Implemented PMoA 6: Decrease in kappa
- Adjusted default params to be compatible with log-logistic responses
