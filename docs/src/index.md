
# AmphiDEB.jl

## Installation

AmphiDEB depends on `EcotoxSystems.jl`, which is curently not registered. 
Installation thus requires

```Julia
using Pkg
Pkg.add(url="htts://github.com/simonhansul/EcotoxSystems.jl")
Pkg.add(url="htts://github.com/simonhansul/AmphiDEB")
```

## TKTD interface

The TKTD module is designed so that mixtures with an arbitrary number of chemicals can be simulated, 
and each chemical can act via arbitrary combinations of PMoAs. <br>
For this pupose, TKTD parameters for sublethal effects are stored in matrices, where the rows represent chemicals and the columns represent PMoAs. <br>
This affects the parameters $k_{D_j}$, $e_j$ and $b_j$. <br>



## API 


```@autodocs
Modules = [AmphiDEB]
Order   = [:type, :function]
```