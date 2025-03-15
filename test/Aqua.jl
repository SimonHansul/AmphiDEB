using Pkg; Pkg.activate("test")
using Aqua
using Revise, AmphiDEB

Aqua.test_all(EcotoxSystems)
