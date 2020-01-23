% EUI - Monetary policy and Inequality under Labor Market Frictions and Capital-Skill Complementarity
% Juan Dolado, Gergõ Motyovszki and Evi Pappa, European University Institute
% codes by Gergõ Motyovszki, EUI, 2017

% running the codes requires IRIS toolbox, 2015jan

% setting parameters and calculation of steady state 

    
%% Parameters

function P = DMP_monpol_CSC_steady(scenarios)
    P = struct();
        
    prodfcn   = scenarios.prodfcn;
    P.prodfcn = scenarios.prodfcn;                  % 1 if complementarities, 0 if CES benchmark
    P.zz      = scenarios.zz;                       % 1 if time-varying capital utilization, 0 if constant
    P.nkpc    = 1;                                  % 1 if non-linear pricing decision, 0 if log-linearized nkpc    
    P.market  = scenarios.market;
    
    
    % Free parameters (27)
    
        % Labor market frictions (8)
            P.alpha_H = scenarios.alpha_H;
            P.alpha_L = scenarios.alpha_L;
            
            P.sigma_H = scenarios.sigma_H;
            P.sigma_L = scenarios.sigma_L;
            P.mu_H = 0;     %scenarios.mu_H;
            P.mu_L = 0;     %scenarios.mu_L;
            P.rho_mIH = scenarios.rho_mIH;
            P.rho_mIL = scenarios.rho_mIL;
            
            
            
        % Household (7)
            P.beta      = scenarios.beta;
            P.eta       = scenarios.eta; 
            P.zeta      = scenarios.zeta;
            P.varphi_H  = scenarios.varphi_H;
            P.varphi_L  = scenarios.varphi_L;
            P.varphi_E  = scenarios.varphi_E;
            P.omega     = scenarios.omega; 
            P.delta     = scenarios.delta; 
            P.h         = scenarios.h;
            
            P.xii       = scenarios.xii;
            P.bH_hat    = scenarios.bH_hat;
            P.bL_hat    = scenarios.bL_hat;
            
            
        % Intermediate firms (6+3)
            P.kappa_H   = scenarios.kappa_H; 
            P.kappa_L   = scenarios.kappa_L; 
            P.a_k       = scenarios.a_k; 
            P.alpha     = scenarios.alpha; 
            P.lambda    = scenarios.lambda; 
            P.gamma     = scenarios.gamma; 
            P.varphi    = scenarios.varphi; 
            P.SK        = scenarios.SK; 
            
            P.sigma = scenarios.sigma; 
            P.xi    = scenarios.xi;
            P.nu    = scenarios.nu; 
            
            
        % Retailers (2)
            P.epsilon   = scenarios.epsilon;
            P.khi       = scenarios.khi;
            
            % endogenized (1)
            %P.tau = ;
            
        % Government (4)
            P.Gamma     = scenarios.Gamma; 
            P.zeta_R    = scenarios.zeta_R; 
            P.zeta_pi   = scenarios.zeta_pi; 
            P.zeta_y    = scenarios.zeta_y; 
            P.zeta_x    = scenarios.zeta_x; 
            P.zeta_u    = scenarios.zeta_u; 
            P.zeta_theta= scenarios.zeta_theta; 

            
        % Autoregressive coeffs (4)
            P.rho_a     = scenarios.rho_a; 
            P.rho_g     = scenarios.rho_g; 
            P.rho_R     = scenarios.rho_R; 
            P.rho_psi   = scenarios.rho_psi; 
            P.rho_inv   = scenarios.rho_inv; 
            P.rho_qq    = scenarios.rho_qq; 
            
            P.rho_wH    = scenarios.rho_wH;
            P.rho_wL    = scenarios.rho_wL;
    
    % Targeted steady state values (10)
    
        
        P.partic_H  = scenarios.partic_H; 
        P.partic_L  = scenarios.partic_L; 
        P.unemp_H   = scenarios.unemp_H; 
        P.unemp_L   = scenarios.unemp_L; 
        P.shrout_H  = 1;    %scenarios.shrout_H; 
        P.shrout_L  = 1;    %scenarios.shrout_L; 
        
        P.bw_H = scenarios.bw_H; 
        P.bw_L = scenarios.bw_L; 
        
        P.z = 1;
        P.x = 1;
        
    % Shock sizes
        P.s_psi = scenarios.s_psi;              % 7 percent increase in the discount rate
        P.s_inv = scenarios.s_inv;              % 5 percent increase in the price of investment
        P.s_a   = scenarios.s_a;
        P.s_g   = scenarios.s_g;
        P.s_R   = scenarios.s_R;
        P.s_qq  = scenarios.s_qq;
        
%% Calculating the Steady State
            
    % Labor market SAM (27)
    
        P.U_H  = P.varphi_H * P.partic_H * P.unemp_H ;     
        P.U_L  = P.varphi_L * P.partic_L * P.unemp_L ;     
        P.rho_mOH = 0;  
        P.rho_mOL = 0;  
        P.N_H   = (1/P.unemp_H - 1) * (P.U_H);
        P.N_L   = (1/P.unemp_L - 1) * (P.U_L);      
        
        P.u_OH = 1;     
        P.u_OL = 1;     
        P.m_IH      = P.sigma_H * P.N_H;   
        P.m_IL      = P.sigma_L * P.N_L;   
                        
        P.v_H   = ( P.m_IH / ( P.rho_mIH * P.U_H^(1-P.alpha_H) ) )^(1/P.alpha_H);
        P.v_L   = ( P.m_IL / ( P.rho_mIL * P.U_L^(1-P.alpha_L) ) )^(1/P.alpha_L);
                
        P.gamma_hIH = P.m_IH / P.U_H;
        P.gamma_hOH = 1;    
        P.gamma_hIL = P.m_IL / P.U_L;
        P.gamma_hOL = 1;    
        
        P.mu_ratio  =  P.gamma_hIH / P.gamma_hIL;
        
        P.gamma_fH       = (P.m_IH) / P.v_H;
        P.gamma_fL       = (P.m_IL) / P.v_L;    
        
        P.theta_H   = P.v_H /(P.U_H ) ;
        P.theta_L   = P.v_L /(P.U_L ) ;
        P.theta     = (P.v_H+P.v_L) / (P.U_H + P.U_L );           
        
        P.m_OH = 1;     
        P.m_OL = 1;     
        
        P.u_IH = P.U_H/P.varphi_H;
        P.u_IL = P.U_L/P.varphi_L;
        P.n_H = P.N_H/P.varphi_H;
        P.n_L = P.N_L/P.varphi_L;
        
    % Households (14)
    
        P.l_H   = 1 - P.n_H - P.u_IH;
        P.l_L   = 1 - P.n_L - P.u_IL;
        
        P.partic  = 1 - P.varphi_H*P.l_H - P.varphi_L*P.l_L ;
        P.unemp   = (P.U_H + P.U_L)/(P.partic) ;
        P.shrout  = 1;  
        
        P.nonemp_H = P.l_H + P.u_IH ;
        P.nonemp_L = P.l_L + P.u_IL ;
        P.nonemp = 1 - P.N_H - P.N_L - P.varphi_E;

        
        P.xbar = P.x;
        P.ubar = P.unemp;
        P.thetabar =P.theta;
        
        
        P.betat = P.beta;
        P.r     = (1-P.beta)/P.beta + P.delta;
        P.Delta = 1;
        P.phi   = P.r / P.delta;
        
        P.lambda_uH = 1;
        P.lambda_uL = 1;
        
    % Pricing and bonds (8)
        
       P.Pi     = 1;                
       P.pi     = log(P.Pi);
       P.R      = P.Pi / P.beta;
       P.rr     = P.R / P.Pi;
       P.tau    = 1/P.epsilon;
       P.qq     = 1;
       
       P.bH     = P.bH_hat;
       P.bL     = P.bL_hat;

       
    % Output and wages (13)
    
        P.A     = 1;
        P.F_k   = P.r / P.x;

%% Basline production function

     if prodfcn == 1 || prodfcn == 2
        display('Production function with complementarities')
        
        % calculating effective capital st st
        
            nH = P.N_H;
            nL = P.N_L;
            AA = P.A;
            ak = P.a_k;
            lambd = P.lambda;
            gamm = P.gamma;
            alph = P.alpha;
            varph = P.varphi;
            if prodfcn == 2
                sk = P.varphi*P.SK;
            else
                sk = P.varphi*P.SK^(1-P.varphi);
            end
            Fk = P.F_k;         
           
            % F_k function (to plot):
            Fkk_fun = @(kt)  sk*ak*lambd*AA * ( ak*(lambd * kt^(gamm) + (1-lambd)*nH^(gamm) )^(alph/gamm) + (1-ak)*nL^(alph) )^(varph*1/alph-1) * ( lambd*kt^(gamm) + (1-lambd)*nH^(gamm) )^((alph-gamm)/gamm) * kt^(gamm-1);
            
                        
            % solve with the symbolic toolbox
            syms kt
            Fk_fun = @(kt) Fk == sk*ak*lambd*AA * ( ak*(lambd * kt^(gamm) + (1-lambd)*nH^(gamm) )^(alph/gamm) + (1-ak)*nL^(alph) )^(varph*1/alph-1) * ( lambd*kt^(gamm) + (1-lambd)*nH^(gamm) )^((alph-gamm)/gamm) * kt^(gamm-1);
            display('Solving for stst effective capital...')
            XX = solve( Fk_fun(kt), kt );
            
            
            P.K = eval(XX);
            
            clear nH nL AA ak lambd gamm alph varph sk Fk kt
                
        P.k     = P.K /( P.z * P.varphi_E);
        if prodfcn == 2
            P.y     =                             P.SK * P.A * ( P.a_k * (P.lambda*P.K^(P.gamma) + (1-P.lambda)*P.N_H^(P.gamma) )^(P.alpha/P.gamma) + (1-P.a_k)*P.N_L^(P.alpha) )^(1/P.alpha *P.varphi);
            P.F_nH  = P.varphi*P.a_k*(1-P.lambda)*P.SK * P.A * ( P.a_k * (P.lambda*P.K^(P.gamma) + (1-P.lambda)*P.N_H^(P.gamma) )^(P.alpha/P.gamma) + (1-P.a_k)*P.N_L^(P.alpha) )^(1/P.alpha *P.varphi -1) * ...
                        (P.lambda * P.K^(P.gamma) + (1-P.lambda)*P.N_H^(P.gamma) )^((P.alpha-P.gamma)/P.gamma) * P.N_H^(P.gamma-1);
            P.F_nL  = P.varphi*(1-P.a_k)         *P.SK * P.A * ( P.a_k * (P.lambda*P.K^(P.gamma) + (1-P.lambda)*P.N_H^(P.gamma) )^(P.alpha/P.gamma) + (1-P.a_k)*P.N_L^(P.alpha) )^(1/P.alpha *P.varphi -1) * P.N_L^(P.alpha-1);
        else
            P.y     =                             P.SK^(1-P.varphi) * P.A * ( P.a_k * (P.lambda*P.K^(P.gamma) + (1-P.lambda)*P.N_H^(P.gamma) )^(P.alpha/P.gamma) + (1-P.a_k)*P.N_L^(P.alpha) )^(1/P.alpha *P.varphi);
            P.F_nH  = P.varphi*P.a_k*(1-P.lambda)*P.SK^(1-P.varphi) * P.A * ( P.a_k * (P.lambda*P.K^(P.gamma) + (1-P.lambda)*P.N_H^(P.gamma) )^(P.alpha/P.gamma) + (1-P.a_k)*P.N_L^(P.alpha) )^(1/P.alpha *P.varphi -1) * ...
                        (P.lambda * P.K^(P.gamma) + (1-P.lambda)*P.N_H^(P.gamma) )^((P.alpha-P.gamma)/P.gamma) * P.N_H^(P.gamma-1);
            P.F_nL  = P.varphi*(1-P.a_k)         *P.SK^(1-P.varphi) * P.A * ( P.a_k * (P.lambda*P.K^(P.gamma) + (1-P.lambda)*P.N_H^(P.gamma) )^(P.alpha/P.gamma) + (1-P.a_k)*P.N_L^(P.alpha) )^(1/P.alpha *P.varphi -1) * P.N_L^(P.alpha-1);
        end
            
        P.wstar_H   = P.x * P.F_nH + (1 - P.sigma_H - 1/P.beta )* P.kappa_H / P.gamma_fH;
        P.wstar_L   = P.x * P.F_nL + (1 - P.sigma_L - 1/P.beta )* P.kappa_L / P.gamma_fL;
        
        P.w_H = P.wstar_H;
        P.w_L = P.wstar_L;
        
        P.b_H   = P.bw_H * P.w_H;
        P.b_L   = P.bw_L * P.w_L;
        
    % Market clearing and aggregates (12)
    
        P.G     = P.y * P.Gamma;
        P.T     = P.b_H*(P.U_H ) + P.b_L*(P.U_L ) + P.G + P.tau*P.x*P.y;
        P.t_H   = P.T;
        P.t_L   = P.T;
        P.t_E   = (P.T - P.varphi_H*P.t_H - P.varphi_L*P.t_L)/P.varphi_E;
        
        P.i     = P.delta * P.k;
        P.Inv   = P.varphi_E * P.i;
        
        P.bE    = (0 - P.varphi_H*P.bH - P.varphi_L*P.bL)/P.varphi_E;
        
        P.C     = P.y - P.Inv - P.G - P.kappa_H*P.v_H - P.kappa_L*P.v_L;
        P.c_H   = P.w_H*P.n_H + P.b_H*P.u_IH - P.t_H - P.bH*(1-P.R/P.Pi);
        P.c_L   = P.w_L*P.n_L + P.b_L*P.u_IL - P.t_L - P.bL*(1-P.R/P.Pi);
        P.c_E   = (P.C - P.varphi_H*P.c_H - P.varphi_L*P.c_L)/P.varphi_E;
                  
    
    % Parameters and Lagrange multipliers   (11)
    
        P.lambda_cH  = ((1-P.h)*P.c_H)^(-P.eta);
        P.lambda_cL  = ((1-P.h)*P.c_L)^(-P.eta);
        P.lambda_cE  = ((1-P.h)*P.c_E)^(-P.eta);
        
        P.FF        = ( P.y * ((1-P.h)*P.c_E)^(-P.eta) ) / ( 1 - P.khi*P.beta );
        P.KK        = P.FF;
              
        
        
        % solving for lambdas and Phi^k                   
            
            ccH=P.lambda_cH; ccL=P.lambda_cL; bet=P.beta; zet=P.zeta;
            bH = P.b_H; 
            bL = P.b_L;             
            lH=P.l_H; lL=P.l_L;
            gamIH=P.gamma_hIH; gamIL=P.gamma_hIL;
            sigH=P.sigma_H; sigL=P.sigma_L;
            wH=P.w_H; wL=P.w_L;
            
            syms nH nL PhH PhL
            
            func_uO_H = @(nH,PhH)       nH == (PhH*lH^(-zet) - bH*ccH  )/gamIH;
            func_n_H  = @(nH,PhH)       nH == bet*(ccH*wH + (1-sigH)*nH  - PhH*lH^(-zet)  );
            
            
            func_uO_L = @(nL,PhL)        nL == (PhL*lL^(-zet) - bL*ccL  )/gamIL;                                       
            func_n_L  = @(nL,PhL)        nL == bet*(ccL*wL + (1-sigL)*nL  - PhL*lL^(-zet)  );                 
            
            
            SS = solve(func_uO_H, func_n_H, func_uO_L, func_n_L, nH,nL,PhH,PhL);
            
         P.Phi_H       = eval(SS.PhH);
         P.Phi_L       = eval(SS.PhL);
         P.lambda_nH   = eval(SS.nH);            
         P.lambda_nL   = eval(SS.nL);           
            
            
            clear bet bH bL gamIH gamIL muH muL nH nL PhH PhL zet
            
        % solving for bargaining powers
        
            nH=P.lambda_nH; nL=P.lambda_nL;
            PhH=P.Phi_H; PhL=P.Phi_L; zet=P.zeta;
            FnH=P.F_nH; FnL=P.F_nL;
            kapH=P.kappa_H; kapL=P.kappa_L;
            gamfH=P.gamma_fH; gamfL=P.gamma_fL;
            
            syms vartH vartL
            
            func_w_H = @(vartH) wH == (1-vartH)*( FnH + (1-sigH)*kapH/gamfH ) + vartH/ccH * (PhH*lH^(-zet) - (1-sigH)*nH  );
            varH = solve(func_w_H, vartH );
         P.vartheta_H = eval(varH);
            
            
            func_w_L = @(vartL) wL == (1-vartL)*( FnL + (1-sigL)*kapL/gamfL ) + vartL/ccL * (PhL*lL^(-zet) - (1-sigL)*nL );
            varL = solve(func_w_L, vartL );
         P.vartheta_L = eval(varL);
                        
            clear cc wH wL sigH sigL lH lL nH nL Ph zet FnH FnL kapH kapL gamfH gamfL vartH vartL
            

    % Taxes and other shocks (2)
            
        P.psi   = 1;
        P.v_R   = 1;
        
    % Other moments to check (1)
        P.kap_v_y = (P.kappa_H * P.v_H + P.kappa_L*P.v_L)/P.y;
        P.kap_wH  = P.kappa_H / P.w_H;
        P.kap_wL  = P.kappa_L / P.w_L;
        P.wprem   = P.w_H / P.w_L;
        P.k_y     = P.K / P.y;
        
    % Income shares and relative prices (13)
        P.skilled = (P.lambda*P.K^(P.gamma) + (1-P.lambda)*(P.N_H)^(P.gamma)  )^(1/P.gamma);
        P.ck = (P.r*P.K + P.w_H*P.N_H) / P.skilled;
        P.r_wH = P.r / P.w_H;
        P.ck_wL = P.ck / P.w_L;
        P.r_wL = P.r / P.w_L;
        P.wH_wL = P.wprem;
        P.share_nH_skilled = P.w_H*P.N_H / ( P.r*P.K + P.w_H*P.N_H );
        P.share_k_skilled  = P.r*P.K     / ( P.r*P.K + P.w_H*P.N_H );
        P.share_skilled_y = P.ck*P.skilled / (P.ck*P.skilled + P.w_L*P.N_L ); 
        P.share_nL_y      = P.w_L*P.N_L    / (P.ck*P.skilled + P.w_L*P.N_L );
        P.share_nH_y      = P.w_H*P.N_H    / (P.ck*P.skilled + P.w_L*P.N_L );
        P.share_k_y       = P.r * P.K      / (P.ck*P.skilled + P.w_L*P.N_L );
        P.rel_labor_share = P.w_H*P.N_H / (P.w_L*P.N_L) ;
        
    
    % Steady states as parameters for the Taylor-rule
        P.Rbar  = P.R;
        P.Pibar = P.Pi;
        P.ybar  = P.y;        
        
    % Other individual vars
       P.profit = (1-(1-P.tau)*P.x)*P.y/P.varphi_E;
       P.borrowing_H = P.rr * P.bH;
       P.borrowing_L = P.rr * P.bL;
       P.borrowing_E = P.rr * P.bE;
       P.net_tr_H = P.b_H*P.u_IH - P.t_H;
       P.net_tr_L = P.b_L*P.u_IL - P.t_L;
       P.net_tr_E = P.tau*P.x*P.y/P.varphi_E - P.t_E;
   
    
    % For LateX tables
    P.StSt1 = [P.l_H; P.l_L; P.n_H; P.n_L; P.u_IH; P.u_IL; P.m_IH; P.m_IL; P.theta_H; P.theta_L; P.bE; P.bH; P.bL];
    P.Stst2 = [P.i; P.k; P.c_E; P.c_H; P.c_L; P.betat; P.R; NaN;  P.lambda_cE; P.lambda_cH; P.lambda_cL; P.lambda_nH; P.lambda_nL];
    P.StSt3 = [P.y; P.K; P.r; P.F_k; P.F_nH; P.F_nL; P.w_H; P.w_L; P.v_H; P.v_L; NaN; P.Pi];            
    P.StSt4 = [P.T; P.psi; P.A; P.G; P.v_R];
    P.Tar1  = [P.gamma_fH; P.gamma_fL; P.gamma_hIH; P.gamma_hIL];
    P.Tar2 =  [P.partic_H; P.partic_L; P.unemp_H; P.unemp_L];
    P.Tar3 =  [P.bw_H; P.bw_L; NaN; P.x; P.z;];
    P.Param1 = [P.sigma_H; P.sigma_L; P.rho_mIH; P.rho_mIL; P.alpha_H; P.alpha_L ];
    P.Param2 = [P.beta; P.eta; P.zeta; P.Phi_H; P.Phi_L; P.varphi_E;P.varphi_H;P.varphi_L; P.omega; P.delta; P.phi; P.h; P.xii];
    P.Param3 = [P.kappa_H; P.kappa_L; P.a_k; P.lambda; P.alpha; P.gamma; P.varphi; 1-P.vartheta_H; 1-P.vartheta_L; P.b_H; P.b_L ];
    P.Param4 = [P.epsilon; P.khi; P.tau; P.Gamma; P.zeta_R; P.zeta_pi; P.zeta_y; P.rho_psi; P.rho_a; P.rho_g; P.rho_R; P.rho_inv ];
    P.Ck1    = [P.k_y; P.kap_v_y; P.kap_wH; P.kap_wL];
    P.Ck2    = [P.partic; P.unemp];
    P.Ck3    = [P.wprem;
                P.rel_labor_share;
                P.share_nH_y;
                P.share_nL_y ;
                P.share_k_y  ;
                P.share_skilled_y ; 
                P.share_nH_skilled;
                P.share_k_skilled];
    

    
%% Benchmark production function
    elseif prodfcn == 0
        display('benchmark CES production function')
        
        % calculating effective capital st st
        
            nH = P.N_H;
            nL = P.N_L;
            AA = P.A;
            sig = P.sigma;
            ksi = P.xi;
            nu = P.nu;
            
            Fk = P.F_k; 
            
            % solve with the symbolic toolbox
            syms kt
            Fk_fun = @(kt) Fk == sig*AA * kt^(sig-1) * ( ksi*nH^nu + (1-ksi)*nL^nu  )^((1-sig)/nu);
            display('Solving for stst effective capital...')
            XX = solve( Fk_fun(kt), kt );
            
            P.K = eval(XX);
            
            clear nH nL AA ak lambd gamm alph varph sk Fk kt
    
        P.k     = P.K / (P.z*P.varphi_E);
        P.y     =                        P.A * P.K^P.sigma * (P.xi*P.N_H^P.nu + (1-P.xi)*P.N_L^P.nu )^((1-P.sigma)/P.nu)  ;
        P.F_nH  =  (1-P.sigma)*P.xi     *P.A * P.K^P.sigma * (P.xi*P.N_H^P.nu + (1-P.xi)*P.N_L^P.nu )^((1-P.sigma-P.nu)/P.nu) * P.N_H^(P.nu-1) ;
        P.F_nL  =  (1-P.sigma)*(1-P.xi) *P.A * P.K^P.sigma * (P.xi*P.N_H^P.nu + (1-P.xi)*P.N_L^P.nu )^((1-P.sigma-P.nu)/P.nu) * P.N_L^(P.nu-1) ;
            
        P.wstar_H   = P.x * P.F_nH + (1 - P.sigma_H - 1/P.beta )* P.kappa_H / P.gamma_fH;
        P.wstar_L   = P.x * P.F_nL + (1 - P.sigma_L - 1/P.beta )* P.kappa_L / P.gamma_fL;
        
        P.w_H = P.wstar_H;
        P.w_L = P.wstar_L;
        
        P.b_H   = P.bw_H * P.w_H;
        P.b_L   = P.bw_L * P.w_L;
        
    % Market clearing and aggregates (12)
    
        P.G     = P.y * P.Gamma;
        P.T     = P.b_H*(P.U_H ) + P.b_L*(P.U_L ) + P.G + P.tau*P.x*P.y;
        P.t_H   = P.T;
        P.t_L   = P.T;
        P.t_E   = (P.T - P.varphi_H*P.t_H - P.varphi_L*P.t_L)/P.varphi_E;
        
        P.i     = P.delta * P.k;
        P.Inv   = P.varphi_E * P.i;
        
        P.bE    = (0 - P.varphi_H*P.bH - P.varphi_L*P.bL)/P.varphi_E;
        
        P.C     = P.y - P.Inv - P.G - P.kappa_H*P.v_H - P.kappa_L*P.v_L;
        P.c_H   = P.w_H*P.n_H + P.b_H*P.u_IH - P.t_H - P.bH*(1-P.R/P.Pi);
        P.c_L   = P.w_L*P.n_L + P.b_L*P.u_IL - P.t_L - P.bL*(1-P.R/P.Pi);
        P.c_E   = (P.C - P.varphi_H*P.c_H - P.varphi_L*P.c_L)/P.varphi_E;
                  
    
    % Parameters and Lagrange multipliers   (11)
    
        P.lambda_cH  = ((1-P.h)*P.c_H)^(-P.eta);
        P.lambda_cL  = ((1-P.h)*P.c_L)^(-P.eta);
        P.lambda_cE  = ((1-P.h)*P.c_E)^(-P.eta);
        
        P.FF        = ( P.y * ((1-P.h)*P.c_E)^(-P.eta) ) / ( 1 - P.khi*P.beta );
        P.KK        = P.FF;
        
        
        
        
        % solving for lambdas and Phi^k                   
            
            ccH=P.lambda_cH; ccL=P.lambda_cL; bet=P.beta; zet=P.zeta;
            bH = P.b_H; 
            bL = P.b_L; 
            lH=P.l_H; lL=P.l_L;
            gamIH=P.gamma_hIH; gamIL=P.gamma_hIL;
            sigH=P.sigma_H; sigL=P.sigma_L;
            wH=P.w_H; wL=P.w_L;
            
            syms nH nL PhH PhL
            
            func_uO_H = @(nH,PhH)       nH == (PhH*lH^(-zet) - bH*ccH  )/gamIH;
            func_n_H  = @(nH,PhH)       nH == bet*(ccH*wH + (1-sigH)*nH  - PhH*lH^(-zet)  );
                        
            func_uO_L = @(nL,PhL)        nL == (PhL*lL^(-zet) - bL*ccL  )/gamIL;                                       
            func_n_L  = @(nL,PhL)        nL == bet*(ccL*wL + (1-sigL)*nL  - PhL*lL^(-zet)  );                 
                        
            SS = solve(func_uO_H, func_n_H, func_uO_L, func_n_L, nH,nL,PhH,PhL);
            
         P.Phi_H       = eval(SS.PhH);
         P.Phi_L       = eval(SS.PhL);
         P.lambda_nH   = eval(SS.nH);            
         P.lambda_nL   = eval(SS.nL);           
            
            
            clear bet bH bL gamIH gamIL muH muL nH nL PhH PhL zet
            
        % solving for bargaining powers
        
            nH=P.lambda_nH; nL=P.lambda_nL;
            PhH=P.Phi_H; PhL=P.Phi_L; zet=P.zeta;
            FnH=P.F_nH; FnL=P.F_nL;
            kapH=P.kappa_H; kapL=P.kappa_L;
            gamfH=P.gamma_fH; gamfL=P.gamma_fL;
            
            syms vartH vartL
            
            func_w_H = @(vartH) wH == (1-vartH)*( FnH + (1-sigH)*kapH/gamfH ) + vartH/ccH * (PhH*lH^(-zet) - (1-sigH)*nH  );
            varH = solve(func_w_H, vartH );
         P.vartheta_H = eval(varH);
            
            
            func_w_L = @(vartL) wL == (1-vartL)*( FnL + (1-sigL)*kapL/gamfL ) + vartL/ccL * (PhL*lL^(-zet) - (1-sigL)*nL );
            varL = solve(func_w_L, vartL );
         P.vartheta_L = eval(varL);
                        
            clear cc wH wL sigH sigL lH lL nH nL Ph zet FnH FnL kapH kapL gamfH gamfL vartH vartL
            

    % Taxes and other shocks (2)
            
        P.psi   = 1;
        P.v_R   = 1;
        
    % Other moments to check (1)
        P.kap_v_y = (P.kappa_H * P.v_H + P.kappa_L*P.v_L)/P.y;
        P.kap_wH  = P.kappa_H / P.w_H;
        P.kap_wL  = P.kappa_L / P.w_L;
        P.wprem   = P.w_H / P.w_L;
        P.k_y     = P.k / P.y;
        
    % Income shares and relative prices (13)
        P.skilled = (P.lambda*P.K^(P.gamma) + (1-P.lambda)*(P.N_H)^(P.gamma)  )^(1/P.gamma);
        P.ck = (P.r*P.K + P.w_H*P.N_H) / P.skilled;
        P.r_wH = P.r / P.w_H;
        P.ck_wL = P.ck / P.w_L;
        P.r_wL = P.r / P.w_L;
        P.wH_wL = P.wprem;
        P.share_nH_skilled = P.w_H*P.N_H / ( P.r*P.K + P.w_H*P.N_H );
        P.share_k_skilled  = P.r*P.K     / ( P.r*P.K + P.w_H*P.N_H );
        P.share_skilled_y = P.ck*P.skilled / (P.ck*P.skilled + P.w_L*P.N_L ); 
        P.share_nL_y      = P.w_L*P.N_L    / (P.ck*P.skilled + P.w_L*P.N_L );
        P.share_nH_y      = P.w_H*P.N_H    / (P.ck*P.skilled + P.w_L*P.N_L );
        P.share_k_y       = P.r * P.K      / (P.ck*P.skilled + P.w_L*P.N_L );
        P.rel_labor_share = P.w_H*P.N_H / (P.w_L*P.N_L) ;
        
    
    % Steady states as parameters for the Taylor-rule
        P.Rbar  = P.R;
        P.Pibar = P.Pi;
        P.ybar  = P.y;        
        
        % Other individual vars
       P.profit = (1-(1-P.tau)*P.x)*P.y/P.varphi_E;
       P.borrowing_H = P.rr * P.bH;
       P.borrowing_L = P.rr * P.bL;
       P.borrowing_E = P.rr * P.bE;
       P.net_tr_H = P.b_H*P.u_IH - P.t_H;
       P.net_tr_L = P.b_L*P.u_IL - P.t_L;
       P.net_tr_E = P.tau*P.x*P.y/P.varphi_E - P.t_E;
        
        
        
    else
        
        error('Wrong value for prodfcn! - It must be either 0 or 1')
        
    end
                
                
                
                
                
                
                
                
                
                
                
                
                