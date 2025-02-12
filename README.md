# AmphiDEB.jl: Modelling amphibian and individual populations with Dynamic Energy Budgets


[![CI](https://github.com/SimonHansul/AmphiDEB.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/SimonHansul/AmphiDEB.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/SimonHansul/AmphiDEB.jl/graph/badge.svg?token=BL1CFR86M6)](https://codecov.io/gh/SimonHansul/AmphiDEB.jl)


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