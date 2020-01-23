% EUI - Monetary policy and Inequality under Labor Market Frictions and Capital-Skill Complementarity
% Juan Dolado, Gergõ Motyovszki and Evi Pappa, European University Institute
% codes by Gergõ Motyovszki, EUI, 2018

% running the codes requires IRIS toolbox, 2015jan

%% Housekeeping
    
    close all
    clear all
    clc
    
        % Loading IRIS Toolbox
%         currentfolder = fileparts(which('DMP_figures.m'));
%         irisfolder = [currentfolder '\IRIS_Tbx_20150127'];
%         addpath(genpath(irisfolder));
% 
%         irisstartup
    
    irisrequired('20150127');
    
    
%% Parameters
    
    baseline = struct();
    
        % Dummy variables
            baseline.prodfcn        = 1;  % 1 with complementarities, 0 with benchmark CES
            baseline.zz             = 0;  % 1 with variable capital utilization, 0 without
            baseline.market         = 3;  % 1: financial autarky (not stable...), 2: incomplete markets, 3: full risk sharing through complete markets
        % Labor market SAM
            baseline.alpha_H        = 0.5;
            baseline.alpha_L        = 0.5;
            baseline.sigma_H        = 0.0245; 
            baseline.sigma_L        = 0.0562; 
            baseline.rho_mIH        = 0.72; 
            baseline.rho_mIL        = 0.455; 
        % Household
            baseline.varphi_H       = 0.21;
            baseline.varphi_L       = 0.69;
            baseline.varphi_E       = 0.10;
            baseline.beta           = 0.99;
            baseline.eta            = 2.0;
            baseline.zeta           = 4;
            baseline.omega          = 4;
            baseline.delta          = 0.01; 
            baseline.h              = 0;
            baseline.bH_hat         = 0.0;
            baseline.bL_hat         = 0.0;
            baseline.xii            = 0.001;
        % Intermediate firms
            baseline.kappa_H        = 0.13;
            baseline.kappa_L        = 0.13;
            baseline.a_k            = 0.4273; 
            baseline.alpha          = 1/2.5;
            baseline.lambda         = 0.35; 
            baseline.gamma          = 1/(-2.04);
            baseline.varphi         = 1;
            baseline.SK             = 1;2.5;
        % Benchmark CES prodfcn
            baseline.sigma          = 0.3; 
            baseline.xi             = 0.5; 
            baseline.nu             = 1; 
        % Retailers
            baseline.epsilon        = 6;
            baseline.khi            = 0.8;
        % Monetary policy and Government
            baseline.Gamma          = 0.2;
            baseline.zeta_R         = 0;
            baseline.zeta_pi        = 1.5;
            baseline.zeta_y         = 0;
            baseline.zeta_x         = 0;
            baseline.zeta_u         = 0;
            baseline.zeta_theta     = 0;
        % Autoregressive coeffs
            baseline.rho_a          = 0.85;
            baseline.rho_g          = 0.7;
            baseline.rho_R          = 0.7;
            baseline.rho_psi        = 0.7;
            baseline.rho_inv        = 0.8;
            baseline.rho_qq         = 0.8;
        % Wage rigidity AR1
            baseline.rho_wH         = 0.0;
            baseline.rho_wL         = 0.0;
        % Targeted steady state values
            baseline.partic_H       = 0.69; 
            baseline.partic_L       = 0.66; 
            baseline.unemp_H        = 0.028; 
            baseline.unemp_L        = 0.078;         
            baseline.bw_H           = 0.4;
            baseline.bw_L           = baseline.bw_H * 1.5306;  % to keep unemployment benefit symmetric instead of replacement rates
        % Shock sizes
            baseline.s_psi = log(0.99);              % 1 percent fall in the discount FACTOR (ppl becoming more impatient, pushing up interest RATES)
            baseline.s_inv = log(0.99);              % 1 percent fall in the price of investment
            baseline.s_a   = log(1.01);
            baseline.s_g   = log(1.01);
            baseline.s_R   = log(0.9975);  % 25 bps cut in the QUARTERLY interest rate (100 bps annually)
            baseline.s_qq  = log(0.99);    % favourable cost-push shock
            
            
            
%% Sensitivity analysis -- different parameters and scenarios
    
        for j=1:80
            scenarios(j) = baseline;
        end

        % 1: remanins the baseline scenario

        % 3: symmetric SAM + Cobb-Douglas
            scenarios(3).partic_H = 0.6;
            scenarios(3).partic_L = 0.6;
            scenarios(3).unemp_H   = 0.06;
            scenarios(3).unemp_L   = 0.06;
            scenarios(3).sigma_H  = 0.0452;
            scenarios(3).sigma_L  = 0.0441;
            scenarios(3).rho_mIH  = 0.6233;
            scenarios(3).rho_mIL  = 0.6131;
            scenarios(3).bw_H  = 0.40; 
            scenarios(3).bw_L  = scenarios(3).bw_H * 1.0000;   
            scenarios(3).kappa_H  = 0.13;
            scenarios(3).kappa_L  = 0.13;
            
            scenarios(3).shrout_H  = 0.16;
            scenarios(3).shrout_L  = 0.16;
            scenarios(3).mu_H  = 0.1333;
            scenarios(3).mu_L  = 0.1333;
            
            scenarios(3).prodfcn = 0;           % standard Cobb-Douglas production function
          
        % 4: symmetric SAM + CSC
            scenarios(4) = scenarios(3);
            scenarios(4).prodfcn = 1;                            % CSC production function
        
        % 5: asymmetric SAM + Cobb-Douglas
            scenarios(5).prodfcn = 0;
            
    % Production function parameters
        % 6-7: capital-skill complementarity -- elasticity
            scenarios(6).gamma      = 1/4;
            scenarios(7).gamma      = -3.5; 
        % 8-9: low skill substitutability -- elasticity
            scenarios(8).alpha    = 1/1.2; %1.2 to make it more substitute
            scenarios(9).alpha    = 1/5;   % making it less substitute  
        
    % labor market SAM friction parameters
        % 14-15: matching efficiencies 
            scenarios(14).rho_mIH = 0.3835;
            scenarios(15).rho_mIH = 0.99;            
        % 16-17: separation rates
            scenarios(16).sigma_H = 0.0505;
            scenarios(17).sigma_H = 0.02;        
        % 26-27: wage rigidities
            scenarios(26).rho_wL = 0.8;%0.3;
            scenarios(26).rho_wH = 0.8;%0.3;
            
            scenarios(27).rho_wL = 0.8;%0.5;
            scenarios(27).rho_wH = 0.6;%0.3;
            
        % 28-29: capital adjustment costs
            scenarios(28).omega      = 8;
            scenarios(29).omega      = 0;                       
        
        % 33-36: monetary policy rules vs CSC
            scenarios(33).zeta_pi   = 100;
            scenarios(34).zeta_y    = 0.5;
            
            scenarios(35).prodfcn = 0;
            scenarios(35).zeta_pi = 100;
            
            scenarios(36).prodfcn = 0;
            scenarios(36).zeta_y = 0.5;
        
        % 37-38: zeta (labor supply elasticity)
            scenarios(37).zeta       = 10;
            scenarios(38).zeta       = 100;                
            
        % 42: variable capital utiliziation
            scenarios(42).zz        = 1;        
            
        
    % labor market SAM friction ONE-BY-ONE
        % 45-46: matching efficiencies 
            scenarios(45) = scenarios(3);
            scenarios(45).rho_mIH = baseline.rho_mIH;
            scenarios(45).rho_mIL = baseline.rho_mIL;
            
            scenarios(46) = scenarios(4);
            scenarios(46).rho_mIH = baseline.rho_mIH;
            scenarios(46).rho_mIL = baseline.rho_mIL;
        % 47-48: separation rates
            scenarios(47) = scenarios(3);
            scenarios(47).sigma_H = baseline.sigma_H;
            scenarios(47).sigma_L = baseline.sigma_L;
            
            scenarios(48) = scenarios(4);
            scenarios(48).sigma_H = baseline.sigma_H;
            scenarios(48).sigma_L = baseline.sigma_L;
        % 49-50: relative weights -- basically bargaining power
            scenarios(49) = scenarios(3);
            scenarios(49).varphi_H = 0.10;%0.36; %to make bargaining symmetric;
            scenarios(49).varphi_L = 0.80;%0.54; %to make bargaining symmetric;
            
            scenarios(50) = scenarios(4);
            scenarios(50).varphi_H = 0.10;
            scenarios(50).varphi_L = 0.80;
            
        % 57-59: variable capital utilization switched on
            scenarios(57) = scenarios(3);  % sym SAM + Cobb Douglas
            scenarios(57).zz = 1;
            
            scenarios(58) = scenarios(4);  % sym SAM + CSC
            scenarios(58).zz = 1;
            
            scenarios(59) = scenarios(5);  % asym SAM + Cobb Douglas
            scenarios(59).zz = 1;
            
            %scenarios(42) is zz = 1 for the baseline asym SAM + CSC
            
        % 77-80) symmetric PHI
            scenarios(77).partic_H = 0.69;
            scenarios(77).partic_L = 0.749;
            
            scenarios(78) = scenarios(3);
            
            scenarios(79) = scenarios(4);
            
            scenarios(80) = scenarios(5);            
            scenarios(80).partic_H = 0.69;
            scenarios(80).partic_L = 0.749;
            
           
        % Scenario labels for figures        
        scens = {'baseline' ,'N/A',  'sym SAM + bnchm', 'sym SAM + CSC', 'as. SAM + bnchm',...
               ['\gamma=' num2str(scenarios(6).gamma)]  ['\gamma=' num2str(scenarios(7).gamma)], ['\alpha=' num2str(scenarios(8).alpha)], ['\alpha=' num2str(scenarios(9).alpha)], 'N/A',...
               'N/A','N/A','N/A',[ '\psi^{H}=' num2str(scenarios(14).rho_mIH)], ['\psi^{H}=' num2str(scenarios(15).rho_mIH)],...
               ['\sigma^H=' num2str(scenarios(16).sigma_H)], ['\sigma^H=' num2str(scenarios(17).sigma_H)], 'N/A','N/A','N/A',...
               'N/A', 'N/A','N/A','N/A','N/A',...
               ['\rho_w^H=\rho_w^L=' num2str(scenarios(26).rho_wH)],['\rho_w^H=' num2str(scenarios(27).rho_wH) ', \rho_w^L=' num2str(scenarios(27).rho_wL)], ['\omega=' num2str(scenarios(28).omega)], ['\omega=' num2str(scenarios(29).omega)],'N/A',...
               'N/A', 'N/A', ['\zeta^{\pi}=' num2str(scenarios(33).zeta_pi)],['\zeta^y=' num2str(scenarios(34).zeta_y)],['CD + \zeta^{\pi}=' num2str(scenarios(35).zeta_pi)],...
               ['CD + \zeta^{y}=' num2str(scenarios(36).zeta_y)],['\xi=' num2str(scenarios(37).zeta)], ['\xi=' num2str(scenarios(38).zeta)], 'N/A', 'N/A',...
               'N/A', 'capital util', 'N/A', 'N/A', 'as. \psi + bnchm',...
               'as. \psi + CSC', 'as. \sigma + bnchm', 'as. \sigma + CSC', 'more as. \vartheta + bnchm', 'more as. \vartheta + CSC',...
               'N/A', 'N/A','N/A','N/A','N/A',...
               'N/A', 'caputil, sym SAM + CD','caputil, sym SAM + CSC','caputil, as.SAM + CD','N/A',...
               'N/A', 'N/A','N/A','N/A','N/A',...
               'N/A', 'N/A','N/A','N/A','N/A',...
               'N/A', 'N/A','N/A','N/A','N/A',...
               'N/A', 'sym \Phi^k: as.SAM + CSC','sym SAM + bnchm','sym SAM + CSC','sym \Phi^k: as.SAM + bnchm'};
               
        % Scenario labels for filenames        
        scen  = {'base' ,'NA2', 'sSAMb', 'sSAMC', 'aSAMb',...
                 'gam1', 'gam2', 'alph1', 'alph2', 'NA10',...
                 'NA11','NA12','NA13','rho1','rho2',...
                 'sig1','sig2','NA18','NA19','NA20',...
                 'NA21','NA22','NA23','NA24','NA25',...
                 'rhow1','rhow2','omeg1','omeg2','NA30',...
                 'NA31','NA32','sIT','uIT','sIT_CD',...
                 'uIT_CD','zet1','zet2','NA39','NA40',...
                 'NA41', 'caput','NA43','NA44', 'scen45',...
                 'scen46','scen47','scen48','scen49','scen50',...
                 'NA51','NA52','NA53','NA54','NA55',...
                 'NA56','NA57','NA58','NA59','NA60',...
                 'NA61','NA62','NA63','NA64','NA65',...
                 'NA66','NA67','NA68','NA69','NA70',...
                 'NA71','NA72','NA73','NA74','NA75',...
                 'NA76','scen77','scen78','scen79','scen80'};
             
       % -------------      
       % Choose which sensitivity analysis scenarios to simulate 
            scennum = [1 3:9 14:17 26:29 33:38 42 45:50 57:59 77:80  ];       % these are the scenarios needed to replicate figures in the paper
       % ------------
             

%% Solving the model - with IRIS toolbox

    for pp=scennum
        
        disp(' ')
        disp(['Running scenario ' num2str(pp) ': '  scens{pp}])
        disp(' ')
         
        % Parameters and steady state
        
            P = DMP_monpol_CSC_steady(scenarios(pp));         % this function sets parameters, targets and solves for the steady state of other variables  
   

        % Model solution

            mod_file = 'DMP_monpol_CSC.model';

            m = model(mod_file,'linear',false,'assign=',P);
       
    
            [flag, discrep, list] = chksstate(m);             % checking whether computed ss satisfies model equations
                if flag
                    disp('Steady state OK');
                else
                    error('Steady state invalid: "%s".\n',list{:});
                end


            [m, npath] = solve(m);                           % solving the model by linearizing it around its steady state (generalized Schur decomposition)
                disp('Solved model');    
                disp(m);
                if npath==1
                    disp('1st order accurate solution OK');
                elseif npath == 0
                    disp('No stable solution, all explosive');
                elseif npath == Inf
                    disp('Multiple stable solution (not unique)');
                end


                % NPath
                %
                % * 1 .. Unique stable solution
                % * 0 .. No stable solution (all explosive)
                % * Inf .. Multiple stable solutions
                % * -1 .. NaN in solved matrices
                % * -2 .. NaN in eigenvalues
                % * -3 .. NaN derivatives in system matrices
                % * -4 .. Steady state does not hold  
   
        
            [T,R,K,Z,H,D,U,Omg] = sspace(m);                 % Calcualting state-space matrices describing the dynamic solution - T is the transition matrix -- for details type help sspace
            eigen = eig(m);                                  % Eigenvalues of the transition matrix
    
        % save current scenario for later
        moodel.(scen{pp}) = m;
        stab.(scen{pp}) = npath;
        stst.(scen{pp}) = flag;
        paramm.(scen{pp}) = P;    
  
    
    end


%% Impulse Response Functions

    NN = 50;
    
    logfun = @log;
    
    disp(' ')
    disp('Simulating Impulse Response Functions...')
    disp(' ')
    
    for pp=scennum
        
        disp(['IRFs for scenario ' num2str(pp) ': '  scens{pp}])
            
        m = moodel.(scen{pp});
        
        tolog = get(m,'logNames');      % getting endogenous model variables to be logged

        % Setting shock sizes
        s_psi = paramm.(scen{pp}).s_psi;        
        s_inv = paramm.(scen{pp}).s_inv;      
        s_a   = paramm.(scen{pp}).s_a;
        s_g   = paramm.(scen{pp}).s_g;
        s_R   = paramm.(scen{pp}).s_R;
        s_qq  = paramm.(scen{pp}).s_qq;

        
        % Discount factor shock 
            dd = zerodb(m,1:NN);          
            dd.shock_psi(1) = s_psi;                     % discount factor shock
            ss = simulate(m,dd,1:NN,'deviation=',true);

            IRF = dboverlay(dd,ss);                                          % x_t / x
            irfs.psi.(scen{pp}) = dbfun(logfun, IRF, 'nameList=', tolog);    % log(x_t/x) ~~ (x_t - x) / x percent deviation from steady state
           
        % Investment price shock 
            dd = zerodb(m,1:NN);          
            dd.shock_inv(1) = s_inv;                     % investment price shock
            ss = simulate(m,dd,1:NN,'deviation=',true);

            IRF = dboverlay(dd,ss);       
            irfs.inv.(scen{pp}) = dbfun(logfun, IRF, 'nameList=', tolog);            

        % TFP shock 
            dd = zerodb(m,1:NN);
            dd.shock_a(1) = s_a;                        % positive productivity shock
            ss = simulate(m,dd,1:NN,'deviation=',true);

            IRF = dboverlay(dd,ss);       
            irfs.a.(scen{pp}) = dbfun(logfun, IRF, 'nameList=', tolog); 

        % Government spending shock 
            dd = zerodb(m,1:NN);
            dd.shock_g(1) = s_g;                        % expansionary fiscal policy
            ss = simulate(m,dd,1:NN,'deviation=',true);

            IRF = dboverlay(dd,ss);       
            irfs.g.(scen{pp}) = dbfun(logfun, IRF, 'nameList=', tolog); 

        % Monetary policy shock 
            dd = zerodb(m,1:NN);
            dd.shock_R(1) = s_R;                        % expansionary monetary policy
            ss = simulate(m,dd,1:NN,'deviation=',true);

            IRF = dboverlay(dd,ss);       
            irfs.v.(scen{pp}) = dbfun(logfun, IRF, 'nameList=', tolog); 
            
        % Cost-push shock 
            dd = zerodb(m,1:NN);
            dd.shock_qq(1) = s_qq;                        % adverse cost-push shock
            ss = simulate(m,dd,1:NN,'deviation=',true);

            IRF = dboverlay(dd,ss);       
            irfs.qq.(scen{pp}) = dbfun(logfun, IRF, 'nameList=', tolog); 

    end   

    sokk = {'psi','a','g','v','dem','inv','qq'};             % shock labels for filenames
    sokkn ={'Discount factor shock', 'TFP shock', 'Gov. exp. shock','Monetary shock', 'Negative demand shock','Investment price shock','Cost-push shock'};  % Shock labels for figures 
    % sokkn = get(m,'ecomments');   
 
        
    
%% ===== %%%%%%%%%%%% ========== %%%%%%%%%%%% ===========    
%%%  Wage bargaining equation - decomposition
% =======================================================

    disp(' ')
    disp('Decomposing the wage bargaining equation...')
    disp(' ')

    IRFS = struct();
    ALPHA = struct();
    TERMS = struct();
        
    wage_scennum = [3 4 5 1 ...
                   57 58 59 42 ...
                   78 79 80 77 ...
                   6 7];                % selecting scenarios to run
     

for pp = wage_scennum
        
    disp(['Decomposition for scenario ' num2str(pp) ': '  scens{pp}])
    
    m   = moodel.(scen{pp});
    irf = irfs.v.(scen{pp});        % selecting the monetary policy shock
    
    deriv = get(m,'derivatives');
    wrt = get(m,'wrt');

    % Parameters
        vartheta_H  = 1-m.vartheta_H;
        vartheta_L  = 1-m.vartheta_L;
        sigma_H     = m.sigma_H;
        sigma_L     = m.sigma_L;
        kappa_H     = m.kappa_H;
        kappa_L     = m.kappa_L;
        Phi_H       = m.Phi_H;
        Phi_L       = m.Phi_L;
        xi          = m.zeta;
        eta         = m.eta;
        varkappa_H  = m.b_H;
        varkappa_L  = m.b_L;
        varsigma    = m.alpha_H;
        psi_H       = m.rho_mIH;
        psi_L       = m.rho_mIL;
        beta        = m.beta;
        
    % Targets
        unemp_H     = m.unemp_H;
        unemp_L     = m.unemp_L;
    
    % Steady states
        x           = m.x;
        F_nH        = m.F_nH;
        F_nL        = m.F_nL;
        w_H         = m.w_H;
        w_L         = m.w_L;
        nu_H        = m.gamma_fH;
        nu_L        = m.gamma_fL;
        lambda_cH   = m.lambda_cH;
        lambda_cL   = m.lambda_cL;
        l_H         = m.l_H;
        l_L         = m.l_L;
        lambda_nH   = m.lambda_nH;
        lambda_nL   = m.lambda_nL;
        mu_H        = m.gamma_hIH;
        mu_L        = m.gamma_hIL;
        
    % IRFs
        x_hat           = irf.x;
        F_nH_hat        = irf.F_nH;
        F_nL_hat        = irf.F_nL;
        l_H_hat         = irf.l_H;
        l_L_hat         = irf.l_L;
        w_H_hat         = irf.w_H;
        w_L_hat         = irf.w_L;
        wprem_hat       = irf.wH_wL;
        theta_H_hat     = irf.theta_H;
        theta_L_hat     = irf.theta_L;
        c_H_hat         = irf.c_H;
        c_L_hat         = irf.c_L;
        
    % log-lin coefficients for the wage bargaining equation
        coeff1_x_H       = vartheta_H * x * F_nH / w_H;
        coeff1_x_L       = vartheta_L * x * F_nL / w_L;
        coeff1_Fn_H      = vartheta_H * x * F_nH / w_H;
        coeff1_Fn_L      = vartheta_L * x * F_nL / w_L;
        coeff1_theta_H   = (1-varsigma) * vartheta_H * (1-sigma_H) * kappa_H / (nu_H * w_H) + varsigma * (1-vartheta_H)*(1-sigma_H)*lambda_nH / (lambda_cH * w_H);
        coeff1_theta_L   = (1-varsigma) * vartheta_L * (1-sigma_L) * kappa_L / (nu_L * w_L) + varsigma * (1-vartheta_L)*(1-sigma_L)*lambda_nL / (lambda_cL * w_L);
        coeff1_c_H       = eta * (1-vartheta_H) / w_H  *  ( ( Phi_H * (l_H)^(-xi) - (1-sigma_H)*lambda_nH )/lambda_cH - (1-sigma_H)*varkappa_H/mu_H );
        coeff1_c_L       = eta * (1-vartheta_L) / w_L  *  ( ( Phi_L * (l_L)^(-xi) - (1-sigma_L)*lambda_nL )/lambda_cL - (1-sigma_L)*varkappa_L/mu_L );
        coeff1_l_H       = xi * (1-vartheta_H) * Phi_H * (l_H)^(-xi) / (lambda_cH * w_H) * ( (1-sigma_H)/mu_H - 1 );
        coeff1_l_L       = xi * (1-vartheta_L) * Phi_L * (l_L)^(-xi) / (lambda_cL * w_L) * ( (1-sigma_L)/mu_L - 1 );
     
        
    % Decomposition of the log-linearized wage bargaining equation
    
        term1_x_H        = coeff1_x_H     * x_hat;
        term1_x_L        = coeff1_x_L     * x_hat;
        term1_Fn_H       = coeff1_Fn_H    * F_nH_hat;
        term1_Fn_L       = coeff1_Fn_L    * F_nL_hat;
        term1_theta_H    = coeff1_theta_H * theta_H_hat;
        term1_theta_L    = coeff1_theta_L * theta_L_hat;
        term1_c_H        = coeff1_c_H     * c_H_hat;
        term1_c_L        = coeff1_c_L     * c_L_hat;
        term1_l_H        = coeff1_l_H     * l_H_hat;
        term1_l_L        = coeff1_l_L     * l_L_hat;
        
        wage1_eq_H       = [term1_x_H term1_Fn_H term1_theta_H term1_c_H term1_l_H];
        wage1_eq_L       = [term1_x_L term1_Fn_L term1_theta_L term1_c_L term1_l_L];
        wage1_prem       = [wage1_eq_H-wage1_eq_L];   

        
    % Saving things
    
        plotlist = {'x_hat', 'F_nH_hat', 'theta_H_hat', 'c_H_hat', 'l_H_hat', 'wprem_hat', 'F_nL_hat', 'theta_L_hat', 'c_L_hat', 'l_L_hat', 'w_H_hat', 'w_L_hat' };
        for vv = 1:length(plotlist)
            IRFS.(scen{pp}).(plotlist{vv}) = eval(plotlist{vv});
        end

        plotlist2 = { 'coeff1_x_H', 'coeff1_x_L', 'coeff1_Fn_H', 'coeff1_Fn_L', 'coeff1_theta_H', 'coeff1_theta_L', 'coeff1_c_H', 'coeff1_c_L', 'coeff1_l_H', 'coeff1_l_L'  };
        for vv = 1:length(plotlist2)
            ALPHA.(scen{pp}).(plotlist2{vv}) = eval(plotlist2{vv});
        end

        plotlist3 = {'wage1_eq_H','wage1_eq_L', 'wage1_prem' };
        for vv = 1:length(plotlist3)
            TERMS.(scen{pp}).(plotlist3{vv}) = eval(plotlist3{vv});
        end
    
end

%% SAVE

    save('DMP_results.mat','moodel','paramm','scenarios','stab','stst','scen','scens','sokk','sokkn','scennum');
    save('DMP_irf_results.mat','irfs','IRFS','ALPHA','TERMS');
