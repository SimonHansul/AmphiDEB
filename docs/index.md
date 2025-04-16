# AmphiDEB.jl 



## Simulating chemical exposure

The global parameters `glb` contain an entry `C_W`, which, by default, is a `Vector{Float64}`:

```Julia
p.glb.C_W = [0.]
```

This represents constant exposure concentrations of all chemical stressors in the system. <br>
In this case, there is one stressor, and the exposure concentration is 0. <br>

In the simplest case, we can simulate different exposure scenarios by modifying 
`p.glb.C_W`. <br><br>

The [EcotoxSystems](https://github.com/simonhansul/ecotoxsystems.jl) package provides a helper function `exposure()`, 
which simplifies this process: 

```Julia

using AmphiDEB, EcotoxSystems

# use the default params
p = deepcopy(AmphiDEB.defaultparams) 

# adjust TKTD parameters 

p.spc.KD .= 0.
p.spc.KD[1,1] = 1.
p.spc.E[1,1] = 1. 
p.spc.B[1,1] = 2.

# use exposure() to simulate a list of treatments and collect the result in a single dataframe

sim = exposure(AmphiDEB.ODE_simulator, p, [0., 1., 2., 4.])

```


### Adding stressors

If we want to add more stressors, we need to re-construct the entire parameter object:

```Julia
p = AmphiDEB.ComponentVector(
    glb = AmphiDEB.ComponentVector(
        t_max = 56., 
        N0 = 1., 
        dX_in = [20., 20.], 
        k_V = [0., 0.], 
        V_patch = [1., 1.], 
        T = 293.15, 
        C_W = [0. 0.;],
        pathogen_inoculation_dose = 0., 
        pathogen_inoculation_time = 30., 
        medium_renewals = [0.] 
    ),
    spc = AmphiDEB.ComponentVector(
        Z = Dirac(1.), 
        propagate_zoom = ComponentVector( 
            dI_max_emb = 1/3, 
            dI_max_lrv = 1/3, 
            dI_max_juv = 1/3, 
            X_emb_int = 1., 
            H_j1 = 1., 
            H_p = 1., 
            K_X_lrv = 1/3, 
            K_X_juv = 1/3
        ),
        
        #=
        Physiological baseline (DEB) parameters
        =#

        X_emb_int = 1, 
        K_X_lrv = 1.,  
        K_X_juv = 1., 
        dI_max_emb = 1, 
        dI_max_lrv = 1, 
        dI_max_juv = 1, 
        kappa_emb = 0.8, 
        kappa_juv = 0.8, 
        gamma = 0.5, 
        eta_IA = 0.54, 
        eta_AS_emb = 0.4, 
        eta_AS_juv = 0.4, 
        eta_AR = 0.95, 
        eta_SA = 0.8, 
        k_M_emb = 0.11, 
        k_M_juv = 0.11, 
        k_J_emb = 0.027, 
        k_J_juv = 0.027, 
        H_j1 = 1, 
        H_p = 55., 

        T_A = 8000., 
        T_ref = 293.15, 
        b_T = 40., 

        #=
        TKTD parameters    
        =#

        h_b = 0.,
        
        KD = [
            0. 0. 0. 0. 0. 0.;
            0. 0. 0. 0. 0. 0.;
            ], 
        B = [
            2. 2. 2. 2. 2. 2.;
            2. 2. 2. 2. 2. 2.;
            ], 
        E = [
            1e10 1e10 1e10 1e10 1e10 1e10;
            1e10 1e10 1e10 1e10 1e10 1e10;
            ], 
        KD_h = [
            0.;
            0.;
            ], 
        E_h = [
            1e10;
            1e10;
            ], 
        B_h = [
            1.;
            1.;
            ], 
        C_h = [
            1.;
            1.;
            ], 

        S_rel_crit = 0.66, 
        h_S = 0.6, 
        a_max = truncated(Normal(15 * 365, 1.5 * 365), 0, Inf), 
        tau_R = 365., 
        
        #=
        Pathogen dynamics and effect parameters
        =#

        Chi = LogNormal(log(1)+1^2, 1), 
        E_P = [Inf, Inf, Inf, Inf], 
        B_P = [2., 2., 2., 2.], 
    )
)
```

Note that the number of entries in `glb.C_W` matches the number of rows in the TKTD parameter matrices, 
since each row in this matrix represents a stressor. <br>

This re-construction is necessary due to how `ComponentArrays.jl` works. <br>
Luckily, it will usually be sufficient to do this construction once in the boilerplate code for a given project and then continue working with (copies of) the default parameters. <br>


## Simulating time-variable exposure

Time-variable exposure is simulated by replacing the scalar entries in `glb` with an interpolator function:

```Julia
p.glb.C_W = [interp_C_W]
```

`interp_C_W` can be an `AbstractInterpolation` or `Function`. <br>
In any case, it has to be possible to call 

```Julia
C_W = interp_C_W(t)
```

to get the interpolated exposure concentration at time `t`. <br><br>

This gives you complete freedom in specifying the interpolation function. <br>
The easier way to simulate time-variable exposure, however, is to again use `exposure()`. <br>
Previously, we had provided `C_Wmat` as a matrix of concentrations to simulate constant exposures. <br>
Now, we provide `C_Wmat`as a `DataFrame`, listing the exposures for all stressors and scenarios. <br>

This dataframe has a fixed format: 

- A column `scenario`, indicating the exposure scenario as numeric value
- A colulmn `t`, referring to time in the simulation 
- A variable number of columns `C_W_1`, `C_W_2`, ... `C_W_n`, listing time-resolved exposure concentrations for stressors $1, 2, ...n$.

The order of the columns does not matter, but the naming does. <br>
The following table specifies a control, two pulsed single-stressor exposures 
and one pulsed mixture exposure. <br>

| scenario | t | C_W_1 | C_W_2 |
|---|---|---|---|
| 0  | 0 | 0 | 0 |
| 0 | 14 | 0 | 0 |
| 1 | 0 | 0 | 0 |
| 1 | 5 | 0 | 0 |
| 1 | 5 | 2 | 0 |
| 1 | 7 | 0 | 0 |
| 1 | 7 | 1 | 0 |
| 1 | 9 | 0 | 0 |
| 1 | 14 | 0 | 0 |
| 2 | 0 | 0 | 0 |
| 2 | 5 | 0 | 0 |
| 2 | 5 | 0 | 0.6 |
| 2 | 7 | 0 | 0 |
| 2 | 7 | 0 | 0.3 |
| 2 | 9 | 0 | 0 |
| 2 | 14 | 0 | 0 |
| 3 | 0 | 0 | 0 |
| 3 | 5 | 0 | 0 |
| 3 | 5 | 2 | 0.6 |
| 3 | 7 | 0 | 0 |
| 3 | 7 | 1 | 0.3 |
| 3 | 9 | 0 | 0 |
| 3 | 14 | 0 | 0 |

To simulate exposure peaks, we can add a zero and non-zero value for the identical time-point. <br>
Attempting to extrapolate beyond the time points listed in the exposure scenario will result in an error. 