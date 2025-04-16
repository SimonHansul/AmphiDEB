import AmphiDEB: ComponentVector

p = ComponentVector(
    glb = ComponentVector(
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
    spc = ComponentVector(
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
    ), 
    pth = AmphiDEB.pth
)