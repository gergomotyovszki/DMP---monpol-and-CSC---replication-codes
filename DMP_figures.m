%% Creating figure for the paper

% EUI - Jobless recoveries and the goals of monetary policy
% Evi Pappa and Juan Dolado, European University Institute
% codes by Gergõ Motyovszki, EUI, 2017

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
    
    % Create folders (if it does not exist) to save figures to
    [~,~,~]=mkdir('numbers');           % printing numerical results
    [~,~,~]=mkdir('charts\paper');      % IRF figures in the main text
    [~,~,~]=mkdir('charts\loglin');     % wage decomposition charts
    [~,~,~]=mkdir('charts\appendix');   % appendix charts

%% Loading results

    load('DMP_irf_results.mat')
    load('DMP_results.mat')

        
%% Printing numbers to text
    
        coefflist     = {'sigma_H', 'sigma_L', 'rho_mIH', 'rho_mIL', 'alpha_H', 'varphi_H', 'varphi_L', 'varphi_E', 'kappa_H', 'b_H',        'b_L',...
                        'rho_a', 'rho_g', 'rho_R', 'lambda', 'alpha', 'gamma', 'omega', 'delta', 'beta', 'eta', 'zeta', 'epsilon', 'khi', 'Gamma', 'zeta_pi', 'zeta_y',...
                        'Phi_H', 'Phi_L', 'vartheta_H', 'vartheta_L', 'tau', 'a_k',...
                        'partic_H', 'partic_L', 'unemp_H', 'unemp_L', 'x', 'wprem',...
                        'theta_H', 'theta_L', 'mu_ratio'};
        coeffnamelist = {'sigma_H', 'sigma_L', 'psi_H',   'psi_L',   'varsigma','varphi_H', 'varphi_L', 'varphi_E', 'kappa',   'varkappa_H', 'varkappa_L',...
                        'rho_a', 'rho_g', 'rho_R', 'lambda', 'alpha', 'gamma', 'omega', 'delta', 'beta', 'eta', 'xi',   'epsilon', 'chi', 'GammaG', 'zeta_pi', 'zeta_y',...
                        'Phi_H', 'Phi_L', 'vartheta_H', 'vartheta_L', 'tau', 'phi',...
                        'partic_H', 'partic_L', 'unemp_H', 'unemp_L', 'x', 'wprem',...
                        'theta_H', 'theta_L', 'mu_ratio'};
        coeffs_printed = zeros(length(coefflist),1);
        
        pp = 1;     % selecting the baseline scenario
               
        
        for i = 1:length(coefflist)
            if i == 30 || i==31         % display bargaining power of workers (instead of of firms)
                coeffs_printed(i) = 1 - paramm.(scen{pp}).(coefflist{i});
            else
                coeffs_printed(i) = paramm.(scen{pp}).(coefflist{i});
            end
            fileID = fopen(['numbers/' coeffnamelist{i} '_' scen{pp} '.tex'],'w');
            fprintf(fileID,'%4.4f',coeffs_printed(i));
        end
        
            
    
%% Setting colors
     
    %Colors
    OR = [232 139 9 ]/255;
    OR2 = [189 94 0]/255;
    GR = [91 92 94]/255;
    GR2= [61 61 61]/255;
    BL = [4  90  209]/255;
    BL2 =[85 154 252]/255;
    myred = [128 0 0]/255;
    myblue= [0 189 189]/255;
    
    barcolormap = [BL;myred;myblue;GR;BL2;OR;GR2;OR2];
    barcolormap1 = [BL;myred;GR;myblue;OR];
    barcoormap6 = [0    0.4470    0.7410;
                0.3010    0.7450    0.9330;
                0.8500    0.3250    0.0980;
                0.4940    0.1840    0.5560;
                0.4660    0.6740    0.1880;
                0.9290    0.6940    0.1250];
    barcoormap5 = [0    0.4470    0.7410;
                0.3010    0.7450    0.9330;
                0.8500    0.3250    0.0980;
                0.4660    0.6740    0.1880;
                0.9290    0.6940    0.1250];
    
    N = 10;
    scale = 100;
    scale_cp = 100;
    
    
%% Scenario constellations to compare

    scenconsts_names = {'Effects of different frictions', 'Variable capital utilization', 'Symmetric labor preferences \Phi^k', 'Different CSC elasticities',...
                        'Different separation rates', 'Different matching efficiencies', 'Different labor supply elasticities', 'Different cap. adj. costs',...
                        'Different elasticities of substitution', 'Different wage rigidities', 'Capital utilization margin', 'Differrent monetary strategies',...
                        'Differrent monetary strategies without CSC','asym \psi + CSC','asym \sigma + CSC','asym \vartheta (pop.weight) + CSC' };
    scenconsts       = {'baseline', 'caputil', 'symphi', 'gamma', 'sigma','psi','xi','omega','alpha','rhow','capz','monpol','monpolCSC','apsi','asigma','avarphi' };
    
        SC.baseline  = [3 4 5 1];           % kk=1
        SC.caputil   = [57 58 59 42];       % kk=2
        SC.symphi    = [78 79 80 77];       % kk=3
        SC.gamma     = [1 6 7];             % kk=4
        SC.sigma     = [1 16 17];           % kk=5
        SC.psi       = [1 14 15];           % kk=6
        SC.xi        = [1 37 38];           % kk=7
        SC.omega     = [1 28 29];           % kk=8
        SC.alpha     = [1 8 9];             % kk=9
        SC.rhow      = [1 26 27];           % kk=10
        SC.capz      = [1 42];              % kk=11
        SC.monpol    = [1 33 34];           % kk=12
        SC.monpolCSC = [5 35 36];           % kk=13
        SC.apsi      = [3 4 45 46];         % kk=14
        SC.asigma    = [3 4 47 48];         % kk=15
        SC.avarphi   = [ 3 4 49 50];        % kk=16
       
        % Creating legend entries
        for ss=1:length(scenconsts)
            leg.(scenconsts{ss}) = scens(SC.(scenconsts{ss}));
        end
        
        % Overwriting some legend entries
        leg.baseline = {'sym SAM + bnchm','sym SAM + CSC','as. SAM + bnchm','as. SAM + CSC'};
        leg.gamma    = {['baseline: 1 / (1-\gamma)=' num2str( round(1/(1-scenarios(1).gamma),2) )],...
                         ['1 / (1-\gamma)=' num2str(1/(1-scenarios(6).gamma) )],...
                         ['1 / (1-\gamma)=' num2str( round(1/(1-scenarios(7).gamma),2) )]};
        leg.monpol    = {'baseline Taylor','strict IT','output react.'};
        leg.monpolCSC = {'Taylor + CD','strict IT + CD','output. react. + CD'};
    
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ===================================================
% =========== FIGURE GENERATION =====================
% ===================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    % Figures of type "f" were created using Matlab 2015b version in the manuscript
    % Figures of type "w" were created using Matlab 2018b (except for w1 for which a commented out section is inserted for these purposes)


%% Figure type f1 - Relative labor share response

    jj = 4;     % monetary policy shock
    pp = 1;     % baseline scenario

    f1=figure();
    
        subplot(1,2,1)
            hh=plot(0:N, scale*[ irfs.(sokk{jj}).(scen{pp}).share_nH_y(0:N)...
                           irfs.(sokk{jj}).(scen{pp}).share_nL_y(0:N)...
                           irfs.(sokk{jj}).(scen{pp}).wH_wL(0:N)], 'linewidth',2.5);
            title('Labor shares of income and skill premium');
            legend('H share','L share','w^H/w^L');
            grid on;
            zeroline;
        
        w_H = resize(irfs.(sokk{jj}).(scen{pp}).w_H, 0:N);
        w_L = resize(irfs.(sokk{jj}).(scen{pp}).w_L, 0:N);
        n_H = resize(irfs.(sokk{jj}).(scen{pp}).n_H, 0:N);
        n_L = resize(irfs.(sokk{jj}).(scen{pp}).n_L, 0:N);     
        
        subplot(1,2,2)
            barcon(0:N, scale*[n_H w_H -n_L -w_L],'colorMap=',barcolormap1);
            legend('N^H','w^H','-N^L','-w^L');
            hold on;
            hh=plot(0:N, scale*[irfs.(sokk{jj}).(scen{pp}).rel_labor_share(0:N) irfs.(sokk{jj}).(scen{pp}).rel_labor_share(0:N)]);
            set(hh(1),'LineWidth',3, 'Color',[1 1 1]);
            set(hh(2),'LineWidth',1.5, 'Color',[0 0 0]);
            title('Relative income share of H vs L');
            grid on;           
            
        filename = ['charts\paper\figure1_' sokk{jj} '_' scen{pp} ];    
        grfun.ftitle(['Income share IRFs to a ' sokkn{jj} ' - ' scens{pp} ]);
            f1=figureFullScreen(f1);
            print(filename,'-dpng'); saveas(f1,filename,'epsc');

    
%% Figure types f2 and f3 - IRFs of main aggregate and relative variables
    
    plotlist2 = {'y', 'i', 'Pi', 'nonemp' };
    namelist2 = {'Output', 'Investment', 'Inflation', 'Non-employment rate' };
    
    plotlist3 = {'nonemp', 'partic', 'wH_wL','rel_labor_share'};
    namelist3 = {'H vs L non-employment rate','H vs L participation rate', 'wage skill premium w^H/w^L','relative income share of H vs L'};
    
    for jj = [1 2 3 4 6 7]                  % shocks: discont factor, TFP, government spending, monetary policy, investment price, cost-push
        
        if jj == 4
            kk_ = [1 3];                    % scenarios to compare: different frictions of asym SAM and CSC under baseline and symmetric Phi calibration
        else
            kk_ = 1;                        % scenarios to compare: different frictions of asym SAM and CSC
        end
        
        if jj == 1
            scale = 100/4;                  % for discount factor shock annualized rate adjustment
        else
            scale = 100;
        end
        
        for kk = kk_                              
        
            scennum = SC.(scenconsts{kk});         
            legef   = leg.(scenconsts{kk});   
        
    f2=figure();    

            for vv = 1:length(plotlist2)
               subplot(2,2,vv);
               ii = 1;
                for pp = scennum
                    hh(ii)=plot(0:N, scale*[irfs.(sokk{jj}).(scen{pp}).(plotlist2{vv})(0:N)], 'linewidth',2.3);
                    hold on;
                    ii=ii+1;
                end
                set(hh(1),'LineStyle','-','Marker','o','MarkerSize',8);
                set(hh(2),'LineStyle','-','Marker','.');
                set(hh(3),'LineStyle','-.','Marker','x','MarkerSize',10);
                set(hh(4),'LineStyle','-.','Marker','.');
                title(namelist2(vv));
                if vv == 2
                    legend(legef,'location','SouthEast');
                else
                end
                grid on;
                zeroline;
            end
            
            grfun.ftitle(['Effects of different frictions - ' sokkn{jj} ]);
            f2=figureFullScreen(f2);
            if kk == 1 && jj == 4
                filename = ['charts\paper\figure2_' sokk{jj} '_' scenconsts{kk}];
                print(filename,'-dpng'); saveas(f2,filename,'epsc');
            else
                filename = ['charts\appendix\figure2_' sokk{jj} '_' scenconsts{kk}];
                print(filename,'-dpng'); saveas(f2,filename,'epsc');
            end
    
     f3=figure();    

            for vv = 1:length(plotlist3)
               subplot(2,2,vv);
               ii = 1;
                for pp = scennum
                    if vv == 1
                        hh(ii)=plot(0:N, scale*[irfs.(sokk{jj}).(scen{pp}).nonemp_H(0:N) - irfs.(sokk{jj}).(scen{pp}).nonemp_L(0:N)], 'linewidth',2);
                        hold on;
                    elseif vv == 2
                        hh(ii)=plot(0:N, scale*[irfs.(sokk{jj}).(scen{pp}).partic_H(0:N) - irfs.(sokk{jj}).(scen{pp}).partic_L(0:N)], 'linewidth',2);
                        hold on;
                    else
                        hh(ii)=plot(0:N, scale*[irfs.(sokk{jj}).(scen{pp}).(plotlist3{vv})(0:N)], 'linewidth',2.3);
                        hold on;
                    end
                    ii=ii+1;
                end
                set(hh(1),'LineStyle','-','Marker','o','MarkerSize',8);
                set(hh(2),'LineStyle','-','Marker','.');
                set(hh(3),'LineStyle','-.','Marker','x','MarkerSize',10);
                set(hh(4),'LineStyle','-.','Marker','.');
                title(namelist3(vv));
                if vv == 4
                    legend(legef,'location','NorthEast');
                else
                end
                grid on;
                zeroline;
            end
            
            grfun.ftitle(['Effects of different frictions - ' sokkn{jj} ]);
            f3=figureFullScreen(f3);
            if kk == 1 && jj == 4
                filename = ['charts\paper\figure3_' sokk{jj} '_' scenconsts{kk} ];
                print(filename,'-dpng'); saveas(f3,filename,'epsc')
            else
                filename = ['charts\appendix\figure3_' sokk{jj} '_' scenconsts{kk} ];
                print(filename,'-dpng'); saveas(f3,filename,'epsc');
            end
            
            
        end 
    end
  

%% Figure type f4 - IRFs of relative labor share with different scenario comparisons

    jj = 4;             % monetary shock
       
    plottitles    = scenconsts_names([1 14 15 16]);             % Adding SAM frictions one-by-one to a symmetric SAM benchmark
    plottitles(1) = {'asym SAM + CSC'};
    
    f4=figure();    
            
            vv = 1;
            for kk = [1 14 15 16]
        
                scennum = SC.(scenconsts{kk});         
                legef   = leg.(scenconsts{kk});                
                 
                subplot(2,2,vv);
                ii = 1;
                    for pp = scennum
                        hh(ii)=plot(0:N, scale*[irfs.(sokk{jj}).(scen{pp}).rel_labor_share(0:N)], 'linewidth',2);
                        hold on;
                        ii=ii+1;
                    end
                set(hh(1),'LineStyle','-','Marker','o','MarkerSize',8);
                set(hh(2),'LineStyle','-','Marker','.');
                set(hh(3),'LineStyle','-.','Marker','x','MarkerSize',10);
                set(hh(4),'LineStyle','-.','Marker','.');
                
                legend(legef,'location','NorthEast');
                
                if vv == 4
                    set(hh(3),'Visible','off')
                    legend([hh(1) hh(2) hh(4)], [legef(1:2) legef(4)],'location','NorthEast');
                else
                end
                title(plottitles(vv));
     
                grid on;
                zeroline;
                vv = vv + 1;
            end           
            
            filename = ['charts\paper\figure4_' sokk{jj}];
            grfun.ftitle(['Relative labor income share -- ' sokkn{jj} ]);
            f4=figureFullScreen(f4);
            print(filename,'-dpng'); saveas(f4,filename,'epsc');


%% Figure type f5 - IRFs of four main variables of interest

    plotlist5 = {'y',      'Pi',        'nonemp',                    'rel_labor_share'};
    namelist5 = {'Output', 'Inflation', 'H vs L non-employment rate','relative income share of H vs L' };
    
    
    for jj = [4 7]                  % shocks: monetary policy, cost-push
        
        if jj == 4
            kk_ = [5 6 7 8 9 10 11 12 13];        % scenarios to compare: sensitivity analysis wrt different structural parameters
        else
            kk_ = [12 13];                        % scenarios to compare: different monetary strategies with or without CSC
        end
        
        for kk = kk_                              
        
            scennum = SC.(scenconsts{kk});         
            legef   = leg.(scenconsts{kk}); 
    
            
            f5=figure();    

            for vv = 1:length(plotlist5)
               subplot(2,2,vv);
               ii = 1;
                for pp = scennum
                    if vv == 3
                        hh(ii)=plot(0:N, scale*[irfs.(sokk{jj}).(scen{pp}).nonemp_H(0:N) - irfs.(sokk{jj}).(scen{pp}).nonemp_L(0:N)], 'linewidth',2.3);
                        hold on;
                    else
                        hh(ii)=plot(0:N, scale*[irfs.(sokk{jj}).(scen{pp}).(plotlist5{vv})(0:N)], 'linewidth',2.3);
                        hold on;
                    end
                    ii=ii+1;
                end
                
                if kk == 11
                    set(hh(1),'LineStyle','-','Linewidth',2.5);
                    set(hh(2),'LineStyle','--','Linewidth',2.5);
                elseif kk == 12 || kk == 13
                    set(hh(1),'LineStyle','-','Marker','o','MarkerSize',8);
                    set(hh(2),'LineStyle','-','Marker','.');
                    set(hh(3),'LineStyle','-','Marker','x','MarkerSize',10);
                else
                    set(hh(1),'LineStyle','-','Linewidth',2.5);
                    set(hh(2),'LineStyle','--','Linewidth',2.5);
                    set(hh(3),'LineStyle','-.');
                end
                
                title(namelist5(vv));
                if vv == 2
                    legend(legef,'location','NorthEast');
                else
                end
                grid on;
                zeroline;
            end
            
            filename = ['charts\appendix\figure5_' sokk{jj} '_' scenconsts{kk}];
            grfun.ftitle([ scenconsts_names{kk} ' -- ' sokkn{jj} ]);
            f5=figureFullScreen(f5);
            print(filename,'-dpng'); saveas(f5,filename,'epsc');

        end
    end
    
    return
    
%% Figure types w1 and w2 - Wage equation decomposition

    scennum = [1 3 4 5 42 77]; 
    
    for pp = scennum         
        
         % Figure type w1
         w1=figure();
         
% here for Matlab2015b version    (for better fontsizes)      
%{         
         subplot(2,2,1)
            bb=barcon(0:N, scale*TERMS.(scen{pp}).wage1_eq_H,'colorMap=',barcoormap5);
            hold on;
            hh=plot(0:N, scale*[IRFS.(scen{pp}).w_H_hat IRFS.(scen{pp}).w_H_hat]);
            set(hh(1),'LineWidth',3, 'Color',[1 1 1]);
            set(hh(2),'LineWidth',1.5, 'Color',[0 0 0]);
            title('w^H');
            grid on; 
            
        subplot(2,2,2)
            bb=barcon(0:N, scale*TERMS.(scen{pp}).wage1_eq_L,'colorMap=',barcoormap5);
            legend('x','Fn','\theta','c','l');
            hold on;
            hh=plot(0:N, scale*[IRFS.(scen{pp}).w_L_hat IRFS.(scen{pp}).w_L_hat]);
            set(hh(1),'LineWidth',3, 'Color',[1 1 1],'HandleVisibility','off');
            set(hh(2),'LineWidth',1.5, 'Color',[0 0 0],'HandleVisibility','off');
            title('w^L');
            grid on;
            
        subplot(2,2,3)
            bb=barcon(0:N, scale*TERMS.(scen{pp}).wage1_prem,'colorMap=',barcoormap5);
            hold on;
            hh=plot(0:N, scale*[IRFS.(scen{pp}).wprem_hat IRFS.(scen{pp}).wprem_hat]);
            set(hh(1),'LineWidth',3, 'Color',[1 1 1],'HandleVisibility','off');
            set(hh(2),'LineWidth',1.5, 'Color',[0 0 0],'HandleVisibility','off');
            title('w^H - w^L');
            grid on; hold off;
        
        grfun.ftitle(['Decomposition of wage bargaining -- ' scens{pp} ]);
%}        
% Matlab2015b ends

            subplot(2,2,1)
                bb=barcon(0:N, scale*TERMS.(scen{pp}).wage1_eq_H,'colorMap=',barcoormap5);
                hold on;
                hh=plot(0:N, scale*[IRFS.(scen{pp}).w_H_hat IRFS.(scen{pp}).w_H_hat]);
                set(hh(1),'LineWidth',6, 'Color',[1 1 1]);
                set(hh(2),'LineWidth',3, 'Color',[0 0 0]);
                title('\fontsize{20} w^H');
                grid on; ax=gca; ax.FontSize = 20;

            subplot(2,2,2)
                bb=barcon(0:N, scale*TERMS.(scen{pp}).wage1_eq_L,'colorMap=',barcoormap5);
                legend('x','Fn','\theta','c','l');
                hold on;
                hh=plot(0:N, scale*[IRFS.(scen{pp}).w_L_hat IRFS.(scen{pp}).w_L_hat]);
                set(hh(1),'LineWidth',6, 'Color',[1 1 1],'HandleVisibility','off');
                set(hh(2),'LineWidth',3, 'Color',[0 0 0],'HandleVisibility','off');
                title('\fontsize{20} w^L');
                grid on; ax=gca; ax.FontSize = 20;

            subplot(2,2,3)
                bb=barcon(0:N, scale*TERMS.(scen{pp}).wage1_prem,'colorMap=',barcoormap5);
                hold on;
                hh=plot(0:N, scale*[IRFS.(scen{pp}).wprem_hat IRFS.(scen{pp}).wprem_hat]);
                set(hh(1),'LineWidth',6, 'Color',[1 1 1],'HandleVisibility','off');
                set(hh(2),'LineWidth',3, 'Color',[0 0 0],'HandleVisibility','off');
                title('\fontsize{20} w^H - w^L');
                grid on; hold off; ax=gca; ax.FontSize = 20;

            grfun.ftitle(['\fontsize{20} Decomposition of wage bargaining -- ' scens{pp} ]);            
            w1=figureFullScreen(w1);
            filename = ['charts\loglin\w1_' scen{pp}];
            print(filename,'-dpng'); saveas(w1,filename,'epsc');
            if pp == 1
                filename = ['charts\paper\w1_' scen{pp}];
                print(filename,'-dpng'); saveas(w1,filename,'epsc');
            else
            end
                    
         % Figure type w2
         w2=figure();

            cc = categorical({'H','L'});
            subplot(2,5,1)
                bb=plot(0:N, scale*IRFS.(scen{pp}).x_hat,'LineWidth',3);
                title('\fontsize{20} x');
                grid on; ax=gca; ax.FontSize = 20;
                zeroline;

            subplot(2,5,6)
                bb=bar(cc,[ ALPHA.(scen{pp}).coeff1_x_H;  ALPHA.(scen{pp}).coeff1_x_L ] );
                bb.FaceColor = 'flat'; bb.CData(2,:) = [0.8500 0.3250 0.0980];
                title('\fontsize{20} \alpha_x');
                grid on; ax=gca; ax.FontSize = 20;

            subplot(2,5,2)
                bb=plot(0:N, scale*[IRFS.(scen{pp}).F_nH_hat IRFS.(scen{pp}).F_nL_hat],'LineWidth',3);
                %legend('H','L');
                title('\fontsize{20} Fn');
                grid on; ax=gca; ax.FontSize = 20;
                zeroline;

            subplot(2,5,7)
                bb=bar(cc,[ ALPHA.(scen{pp}).coeff1_Fn_H;  ALPHA.(scen{pp}).coeff1_Fn_L ] );
                bb.FaceColor = 'flat'; bb.CData(2,:) = [0.8500 0.3250 0.0980];
                title('\fontsize{20} \alpha_{Fn}');
                grid on; ax=gca; ax.FontSize = 20;

            subplot(2,5,3)
                bb=plot(0:N, scale*[IRFS.(scen{pp}).theta_H_hat IRFS.(scen{pp}).theta_L_hat],'LineWidth',3);
                legend('H','L');
                title('\fontsize{20} \theta');
                grid on; ax=gca; ax.FontSize = 20;
                zeroline;

            subplot(2,5,8)
                bb=bar(cc,[ ALPHA.(scen{pp}).coeff1_theta_H;  ALPHA.(scen{pp}).coeff1_theta_L ] );
                bb.FaceColor = 'flat'; bb.CData(2,:) = [0.8500 0.3250 0.0980];
                title('\fontsize{20} \alpha_{\theta}');
                grid on; ax=gca; ax.FontSize = 20;

            subplot(2,5,4)
                bb=plot(0:N, scale*[IRFS.(scen{pp}).c_H_hat IRFS.(scen{pp}).c_L_hat],'LineWidth',3);
                title('\fontsize{20} c');
                grid on; ax=gca; ax.FontSize = 20;
                zeroline;

            subplot(2,5,9)
                bb=bar(cc,[ ALPHA.(scen{pp}).coeff1_c_H;  ALPHA.(scen{pp}).coeff1_c_L ] );
                bb.FaceColor = 'flat'; bb.CData(2,:) = [0.8500 0.3250 0.0980];
                title('\fontsize{20} \alpha_{c}');
                grid on; ax=gca; ax.FontSize = 20;

            subplot(2,5,5)
                bb=plot(0:N, scale*[IRFS.(scen{pp}).l_H_hat IRFS.(scen{pp}).l_L_hat],'LineWidth',3);
                title('\fontsize{20} l');
                grid on; ax=gca; ax.FontSize = 20;
                zeroline;

            subplot(2,5,10)
                bb=bar(cc,[ ALPHA.(scen{pp}).coeff1_l_H;  ALPHA.(scen{pp}).coeff1_l_L ] );
                bb.FaceColor = 'flat'; bb.CData(2,:) = [0.8500 0.3250 0.0980];
                title('\fontsize{20} \alpha_{l}');
                grid on; ax=gca; ax.FontSize = 20;

            filename = ['charts\loglin\w2_' scen{pp} ];
            grfun.ftitle(['\fontsize{20} IRFs and wage prem coefficients -- ' scens{pp} ]);
            w2=figureFullScreen(w2);
            print(filename,'-dpng'); saveas(w2,filename,'epsc');

    end

    
%% Figure types w3, w4 and w5 -- Comparing wage decompositions across scenarios

    plotlist_w3 = { 'coeff1_x_H', 'coeff1_x_L', 'coeff1_Fn_H', 'coeff1_Fn_L', 'coeff1_theta_H', 'coeff1_theta_L', 'coeff1_c_H', 'coeff1_c_L', 'coeff1_l_H', 'coeff1_l_L'  };
    
    plotlist_w4 = {'x_hat', 'F_nH_hat', 'theta_H_hat', 'c_H_hat', 'l_H_hat', 'wprem_hat', 'F_nL_hat', 'theta_L_hat', 'c_L_hat', 'l_L_hat' };
    namelist_w4 = {'x', 'Fn^H', '\theta^H', 'c^H', 'l^H', 'w^H - w^L', 'Fn^L', '\theta^L', 'c^L', 'l^L' };         
       
        
    for kk = [1 2 3 4]
        
        scennum = SC.(scenconsts{kk});
        legef   = leg.(scenconsts{kk}); 
        
        
        % Figure type w5
        w5 = figure();
            
            lll = 1;
            for ll = [1 6]  % selecting variables from plotlist_w4: marginal costs (x) and wage premium (wprem)
                
            subplot(2,2,2*(lll-0.5))
                plotmatrix = tseries(0:N,zeros(N+1,length(scennum)));
                ii=1;
                for pp=scennum
                    plotmatrix(:,ii) = [IRFS.(scen{pp}).(plotlist_w4{ll})];
                    ii=ii+1;
                end

                bb=plot(0:N, scale*plotmatrix,'LineWidth',3);
                if kk == 4
                    set(bb(1),'LineStyle','-','Linewidth',3);
                    set(bb(2),'LineStyle','--','Linewidth',3);
                    set(bb(3),'LineStyle','-.');
                else
                    set(bb(1),'LineStyle','-','Marker','o','MarkerSize',12);
                    set(bb(2),'LineStyle','-','Marker','.');
                    set(bb(3),'LineStyle','-.','Marker','x','MarkerSize',16);
                    set(bb(4),'LineStyle','-.','Marker','.');
                end
                title(['\fontsize{20} ' namelist_w4{ll}]);
                grid on; ax=gca; ax.FontSize = 20;
                zeroline;
                if lll == 1
                    legend( legef )
                else
                end
                
            lll = lll+1;    
            end
                
           cc = categorical({'H','L'});
            subplot(2,2,2)
                barmatrix = zeros(2,length(scennum));
                ii=1;
                for pp=scennum
                    barmatrix(:,ii) = [ALPHA.(scen{pp}).(plotlist_w3{1}); ALPHA.(scen{pp}).(plotlist_w3{2}) ];
                    ii=ii+1;
                end
                bb=bar(cc, barmatrix );
                title('\fontsize{20} \alpha_{x}');
                grid on; ax=gca; ax.FontSize = 20;

           cc = categorical({'H-L',' '});
           subplot(2,2,4)
                barmatrix = zeros(2,length(scennum));
                ii=1;
                for pp=scennum
                    barmatrix(1,ii) = [ALPHA.(scen{pp}).(plotlist_w3{1}) - ALPHA.(scen{pp}).(plotlist_w3{2}) ];
                    ii=ii+1;
                end
                bb=bar(cc, barmatrix );
                title('\fontsize{20} \alpha_{x}^H - \alpha_{x}^L');
                grid on; ax=gca; ax.FontSize = 20;

            filename = ['charts\loglin\w5_' scenconsts{kk}];  % filename = ['charts\loglin\w5_' scenconsts{kk} '_fix'];
            grfun.ftitle('\fontsize{20} IRFs and wage coeffs (monetary shock)');
            w5=figureFullScreen(w5);
            print(filename,'-dpng'); saveas(w5,filename,'epsc'); % saveas(gcf,filename,'epsc');
            if kk == 1 || kk == 4
                filename = ['charts\paper\w5_' scenconsts{kk}];
                print(filename,'-dpng'); saveas(w5,filename,'epsc');
            else
            end
            
      
        if kk == 1 || kk == 2        % Do w3 and w4 type figures only for kk = 1 or kk = 2
            
        % Figure type w3
        w3 = figure();
        
            cc = categorical({'H','L'});
            subplot(2,3,1)
                barmatrix = zeros(2,length(scennum));
                ii=1;
                for pp=scennum
                    barmatrix(:,ii) = [ALPHA.(scen{pp}).(plotlist_w3{1}); ALPHA.(scen{pp}).(plotlist_w3{2}) ];
                    ii=ii+1;
                end
                bb=bar(cc, barmatrix );
                title('\fontsize{20} \alpha_{x}');
                grid on; ax=gca; ax.FontSize = 20;

            subplot(2,3,2)
                barmatrix = zeros(2,length(scennum));
                ii=1;
                for pp=scennum
                    barmatrix(:,ii) = [ALPHA.(scen{pp}).(plotlist_w3{3}); ALPHA.(scen{pp}).(plotlist_w3{4}) ];
                    ii=ii+1;
                end
                bb=bar(cc, barmatrix );

                legend( legef )            

                title('\fontsize{20} \alpha_{Fn} = \alpha_x');
                grid on; ax=gca; ax.FontSize = 20;

            subplot(2,3,3)
                barmatrix = zeros(2,length(scennum));
                ii=1;
                for pp=scennum
                    barmatrix(:,ii) = [ALPHA.(scen{pp}).(plotlist_w3{5}); ALPHA.(scen{pp}).(plotlist_w3{6}) ];
                    ii=ii+1;
                end
                bb=bar(cc, barmatrix );
                title('\fontsize{20} \alpha_{\theta}');
                grid on; ax=gca; ax.FontSize = 20;

            subplot(2,3,5)
                barmatrix = zeros(2,length(scennum));
                ii=1;
                for pp=scennum
                    barmatrix(:,ii) = [ALPHA.(scen{pp}).(plotlist_w3{7}); ALPHA.(scen{pp}).(plotlist_w3{8}) ];
                    ii=ii+1;
                end
                bb=bar(cc, barmatrix );
                title('\fontsize{20} \alpha_{c}');
                grid on; ax=gca; ax.FontSize = 20;

            subplot(2,3,6)
                barmatrix = zeros(2,length(scennum));
                ii=1;
                for pp=scennum
                    barmatrix(:,ii) = [ALPHA.(scen{pp}).(plotlist_w3{9}); ALPHA.(scen{pp}).(plotlist_w3{10}) ];
                    ii=ii+1;
                end
                bb=bar(cc, barmatrix );
                title('\fontsize{20} \alpha_{l}');
                grid on; ax=gca; ax.FontSize = 20;

            cc = categorical({'H-L',' '});
            subplot(2,3,4)
                barmatrix = zeros(2,length(scennum));
                ii=1;
                for pp=scennum
                    barmatrix(1,ii) = [ALPHA.(scen{pp}).(plotlist_w3{1}) - ALPHA.(scen{pp}).(plotlist_w3{2}) ];
                    ii=ii+1;
                end
                bb=bar(cc, barmatrix );
                title('\fontsize{20} \alpha_{x}^H - \alpha_{x}^L');
                grid on; ax=gca; ax.FontSize = 20;

            filename = ['charts\loglin\w3_' scenconsts{kk}];
            grfun.ftitle('\fontsize{20} Coefficients of the wage bargaining eq.');
            w3=figureFullScreen(w3);        
            print(filename,'-dpng'); saveas(w3,filename,'epsc');
        
            
        % Figure type w4
         w4 = figure();

            for vv=1:length(plotlist_w4)

                subplot(2,5,vv)
                    plotmatrix = tseries(0:N,zeros(N+1,length(scennum)));
                    ii=1;
                    for pp=scennum
                        plotmatrix(:,ii) = [IRFS.(scen{pp}).(plotlist_w4{vv})];
                        ii=ii+1;
                    end
                    bb=plot(0:N, scale*plotmatrix,'LineWidth',3);
                    title(['\fontsize{20} ' namelist_w4{vv}]);
                    grid on; ax=gca; ax.FontSize = 20;
                    zeroline;

            end
            legend( legef );

            filename = ['charts\loglin\w4_' scenconsts{kk}];
            grfun.ftitle('\fontsize{20} IRFs to a monetary shock');
            w4=figureFullScreen(w4);
            print(filename,'-dpng'); saveas(w4,filename,'epsc');
  
        
        else
        end
        
 
    end   
   