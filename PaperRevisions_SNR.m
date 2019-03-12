%%% Performs MC simulations for (new) Supplementary Figure 1. %%%

close all; clear all;

T1_F = 0.45; T1_S = 1.4; T2_F = 0.015; T2_S = 0.09; M0_F = 0.15; k_FS = 8; Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;
TR_SPGR = 6.5e-3; TR_SSFP = 6.5e-3; FA_SPGR = deg2rad([2 4 6 8 10 12 14 16 18 20]); FA_SSFP180 = deg2rad([2 6 14 22 30 38 46 54 62 70]); FA_SSFP0 = deg2rad([2 6 14 22 30 38 46 54 62 70]);
    
%% Perform stochastic region contraction.

Realisations = 1000; Trials = 40000; Iterations = 30; N = 50; Runs = 1; Params = 6;
 
% Ground-truth signals for mcDESPOT.
SPGR_Data = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'M0_F',M0_F,'k_FS',k_FS);
SSFP_Data_0 = SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);
SSFP_Data_180 = SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);

% Concatenate SSFP signals.
SSFP_Data = [SSFP_Data_0 ; SSFP_Data_180];

% Define SNR and define wrt mean SPGR signal.
SNR = [30,50,100,200];

Solution_Bouhrara = zeros(length(SNR),Realisations,Params);

delete(gcp('nocreate')); c = parcluster('local'); c.NumWorkers = 8; parpool(c, c.NumWorkers);

for ii = 1:length(SNR)
    
    Sigma = mean(SPGR_Data)/SNR(ii);
    
    parfor tt = 1:Realisations
        
        % Add different noise to each element.
        SPGR_Data_Noisy = zeros(length(SPGR_Data),1);
        for mm = 1:length(SPGR_Data)
            SPGR_Data_Noisy(mm) = SPGR_Data(mm) + (normrnd(0,Sigma));
        end
        SSFP_Data_Noisy = zeros(length(SSFP_Data),1);
        for nn = 1:length(SSFP_Data)
            SSFP_Data_Noisy(nn) = SSFP_Data(nn) + (normrnd(0,Sigma));
        end
        % Normalise signals and concatenate.
        SPGR_Data_Norm = SPGR_Data_Noisy./mean(SPGR_Data_Noisy);
        SSFP_Data_Norm = SSFP_Data_Noisy./mean(SSFP_Data_Noisy);
        Data = [SPGR_Data_Norm ; SSFP_Data_Norm];
        
        [Solution_Bouhrara(ii,tt,:), ~, ~] = SRC_SimB(Trials, Iterations, N, Runs, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data);
        
    end
    
end

%%

figure(1)
cm = colormap(lines(4));
subplot(2,3,1);
histogram(Solution_Bouhrara(4,:,1),'BinWidth',0.05,'EdgeColor',cm(1,:),'DisplayStyle','stairs','LineWidth',2); hold on;
histogram(Solution_Bouhrara(3,:,1),'BinWidth',0.05,'EdgeColor',cm(2,:),'DisplayStyle','stairs','LineWidth',2);
histogram(Solution_Bouhrara(2,:,1),'BinWidth',0.05,'EdgeColor',cm(3,:),'DisplayStyle','stairs','LineWidth',2); 
histogram(Solution_Bouhrara(1,:,1),'BinWidth',0.05,'EdgeColor',cm(4,:),'DisplayStyle','stairs','LineWidth',2);
line([T1_S T1_S],[0 400],'LineStyle','--','Color','k','LineWidth',2) 
xlabel('T_{1S} (s)','FontSize',20); ylabel('Count','FontSize',20); 
get(gca, 'XTick'); set(gca, 'FontSize', 20); get(gca, 'YTick'); set(gca, 'FontSize', 20);
grid on; grid minor;
ll = legend('200','100','50','30'); ll.FontSize = 18; title(ll,'SNR','FontSize',22);

subplot(2,3,2);
histogram(Solution_Bouhrara(4,:,2),'BinWidth',0.02,'EdgeColor',cm(1,:),'DisplayStyle','stairs','LineWidth',2); hold on;
histogram(Solution_Bouhrara(3,:,2),'BinWidth',0.02,'EdgeColor',cm(2,:),'DisplayStyle','stairs','LineWidth',2);
histogram(Solution_Bouhrara(2,:,2),'BinWidth',0.02,'EdgeColor',cm(3,:),'DisplayStyle','stairs','LineWidth',2); 
histogram(Solution_Bouhrara(1,:,2),'BinWidth',0.02,'EdgeColor',cm(4,:),'DisplayStyle','stairs','LineWidth',2);
line([T1_F T1_F],[0 300],'LineStyle','--','Color','k','LineWidth',2) 
xlabel('T_{1F} (s)','FontSize',20); ylabel('Count','FontSize',20); 
get(gca, 'XTick'); set(gca, 'FontSize', 20); get(gca, 'YTick'); set(gca, 'FontSize', 20);
grid on; grid minor;

subplot(2,3,3);
histogram(Solution_Bouhrara(4,:,3),'BinWidth',0.005,'EdgeColor',cm(1,:),'DisplayStyle','stairs','LineWidth',2); hold on;
histogram(Solution_Bouhrara(3,:,3),'BinWidth',0.005,'EdgeColor',cm(2,:),'DisplayStyle','stairs','LineWidth',2);
histogram(Solution_Bouhrara(2,:,3),'BinWidth',0.005,'EdgeColor',cm(3,:),'DisplayStyle','stairs','LineWidth',2); 
histogram(Solution_Bouhrara(1,:,3),'BinWidth',0.005,'EdgeColor',cm(4,:),'DisplayStyle','stairs','LineWidth',2);
line([T2_S T2_S],[0 500],'LineStyle','--','Color','k','LineWidth',2) 
xlabel('T_{2S} (s)','FontSize',20); ylabel('Count','FontSize',20); 
get(gca, 'XTick'); set(gca, 'FontSize', 20); get(gca, 'YTick'); set(gca, 'FontSize', 20);
grid on; grid minor;

subplot(2,3,4);
histogram(Solution_Bouhrara(4,:,4),'BinWidth',0.001,'EdgeColor',cm(1,:),'DisplayStyle','stairs','LineWidth',2); hold on;
histogram(Solution_Bouhrara(3,:,4),'BinWidth',0.001,'EdgeColor',cm(2,:),'DisplayStyle','stairs','LineWidth',2);
histogram(Solution_Bouhrara(2,:,4),'BinWidth',0.001,'EdgeColor',cm(3,:),'DisplayStyle','stairs','LineWidth',2); 
histogram(Solution_Bouhrara(1,:,4),'BinWidth',0.001,'EdgeColor',cm(4,:),'DisplayStyle','stairs','LineWidth',2);
line([T2_F T2_F],[0 400],'LineStyle','--','Color','k','LineWidth',2) 
xlabel('T_{2F} (s)','FontSize',20); ylabel('Count','FontSize',20); 
get(gca, 'XTick'); set(gca, 'FontSize', 20); get(gca, 'YTick'); set(gca, 'FontSize', 20);
grid on; grid minor;

subplot(2,3,5);
histogram(Solution_Bouhrara(4,:,5),'BinWidth',0.02,'EdgeColor',cm(1,:),'DisplayStyle','stairs','LineWidth',2); hold on;
histogram(Solution_Bouhrara(3,:,5),'BinWidth',0.02,'EdgeColor',cm(2,:),'DisplayStyle','stairs','LineWidth',2);
histogram(Solution_Bouhrara(2,:,5),'BinWidth',0.02,'EdgeColor',cm(3,:),'DisplayStyle','stairs','LineWidth',2); 
histogram(Solution_Bouhrara(1,:,5),'BinWidth',0.02,'EdgeColor',cm(4,:),'DisplayStyle','stairs','LineWidth',2);
line([M0_F M0_F],[0 500],'LineStyle','--','Color','k','LineWidth',2) 
xlabel('MWF','FontSize',20); ylabel('Count','FontSize',20); 
get(gca, 'XTick'); set(gca, 'FontSize', 20); get(gca, 'YTick'); set(gca, 'FontSize', 20);
grid on; grid minor;

subplot(2,3,6);
histogram(Solution_Bouhrara(4,:,6),'BinWidth',1,'EdgeColor',cm(1,:),'DisplayStyle','stairs','LineWidth',2); hold on;
histogram(Solution_Bouhrara(3,:,6),'BinWidth',1,'EdgeColor',cm(2,:),'DisplayStyle','stairs','LineWidth',2);
histogram(Solution_Bouhrara(2,:,6),'BinWidth',1,'EdgeColor',cm(3,:),'DisplayStyle','stairs','LineWidth',2); 
histogram(Solution_Bouhrara(1,:,6),'BinWidth',1,'EdgeColor',cm(4,:),'DisplayStyle','stairs','LineWidth',2);
line([k_FS k_FS],[0 400],'LineStyle','--','Color','k','LineWidth',2) 
xlabel('k_{FS} (s)','FontSize',20); ylabel('Count','FontSize',20); 
get(gca, 'XTick'); set(gca, 'FontSize', 20); get(gca, 'YTick'); set(gca, 'FontSize', 20);
grid on; grid minor;
