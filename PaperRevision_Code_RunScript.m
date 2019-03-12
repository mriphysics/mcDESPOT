%%% Test of SRC bias. Generates data for (new) Figure 8. %%%

%% Different MWF.

T1_F = 0.45; T1_S = 1.4; T2_F = 0.015; T2_S = 0.09; k_FS = 8; Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;

Params = 6;
TR_SPGR = 6.5e-3; TR_SSFP = 6.5e-3; FA_SPGR = deg2rad([2 4 6 8 10 12 14 16 18 20]); FA_SSFP180 = deg2rad([2 6 14 22 30 38 46 54 62 70]); FA_SSFP0 = deg2rad([2 6 14 22 30 38 46 54 62 70]);

M0_F = linspace(0,0.25,8);

Realisations = 1000; Trials = 40000; Iterations = 30; N = 50; Runs = 1; 

Solution_Bouhrara_MWF = zeros(length(M0_F),Realisations,Params);

delete(gcp('nocreate')); c = parcluster('local'); c.NumWorkers = 8; parpool(c, c.NumWorkers);

for ii = 1:length(M0_F)
    
    % Ground-truth signals for mcDESPOT.
    SPGR_Data = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'M0_F',M0_F(ii),'k_FS',k_FS);
    SSFP_Data_0 = SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F(ii),'k_FS',k_FS);
    SSFP_Data_180 = SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F(ii),'k_FS',k_FS);
    
    % Concatenate SSFP signals.
    SSFP_Data = [SSFP_Data_0 ; SSFP_Data_180];
    
    % Define SNR and define wrt mean SPGR signal.
    SNR = 100; Sigma = mean(SPGR_Data)/SNR;
    
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

        [Solution_Bouhrara_MWF(ii,tt,:), ~, ~] = SRC_SimD(Trials, Iterations, N, Runs, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data);
        
    end
    
end

disp('MWF Done')

%% Different T1.

T1_F = 0.45; T2_F = 0.015; T2_S = 0.09; M0_F = 0.15; k_FS = 8; Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;

Params = 6;
TR_SPGR = 6.5e-3; TR_SSFP = 6.5e-3; FA_SPGR = deg2rad([2 4 6 8 10 12 14 16 18 20]); FA_SSFP180 = deg2rad([2 6 14 22 30 38 46 54 62 70]); FA_SSFP0 = deg2rad([2 6 14 22 30 38 46 54 62 70]);

T1_S = linspace(1.0,2.0,8);

Realisations = 1000; Trials = 40000; Iterations = 30; N = 50; Runs = 1; 

Solution_Bouhrara_T1 = zeros(length(T1_S),Realisations,Params);

delete(gcp('nocreate')); c = parcluster('local'); c.NumWorkers = 8; parpool(c, c.NumWorkers);

for ii = 1:length(T1_S)
    
    % Ground-truth signals for mcDESPOT.
    SPGR_Data = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1_S(ii),'T1_F',T1_F,'M0_F',M0_F,'k_FS',k_FS);
    SSFP_Data_0 = SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1,'T1_S',T1_S(ii),'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);
    SSFP_Data_180 = SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2,'T1_S',T1_S(ii),'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);
    
    % Concatenate SSFP signals.
    SSFP_Data = [SSFP_Data_0 ; SSFP_Data_180];
    
    % Define SNR and define wrt mean SPGR signal.
    SNR = 100; Sigma = mean(SPGR_Data)/SNR;
    
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

        [Solution_Bouhrara_T1(ii,tt,:), ~, ~] = SRC_SimD(Trials, Iterations, N, Runs, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data);
        
    end
    
end

disp('T1 Done')

%% Different T2.

T1_F = 0.45; T1_S = 1.4; T2_F = 0.015; M0_F = 0.15; k_FS = 8; Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;

Params = 6;
TR_SPGR = 6.5e-3; TR_SSFP = 6.5e-3; FA_SPGR = deg2rad([2 4 6 8 10 12 14 16 18 20]); FA_SSFP180 = deg2rad([2 6 14 22 30 38 46 54 62 70]); FA_SSFP0 = deg2rad([2 6 14 22 30 38 46 54 62 70]);

T2_S = linspace(0.04,0.1,8);

Realisations = 1000; Trials = 40000; Iterations = 30; N = 50; Runs = 1; 

Solution_Bouhrara_T2 = zeros(length(T2_S),Realisations,Params);

delete(gcp('nocreate')); c = parcluster('local'); c.NumWorkers = 8; parpool(c, c.NumWorkers);

for ii = 1:length(T2_S)
    
    % Ground-truth signals for mcDESPOT.
    SPGR_Data = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'M0_F',M0_F,'k_FS',k_FS);
    SSFP_Data_0 = SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1,'T1_S',T1_S,'T2_S',T2_S(ii),'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);
    SSFP_Data_180 = SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2,'T1_S',T1_S,'T2_S',T2_S(ii),'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);
    
    % Concatenate SSFP signals.
    SSFP_Data = [SSFP_Data_0 ; SSFP_Data_180];
    
    % Define SNR and define wrt mean SPGR signal.
    SNR = 100; Sigma = mean(SPGR_Data)/SNR;
    
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

        [Solution_Bouhrara_T2(ii,tt,:), ~, ~] = SRC_SimD(Trials, Iterations, N, Runs, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data);
        
    end
    
end

disp('T2 Done')

%%

cm = lines(6);

T1_F = 0.45; T1_S = 1.4; T2_F = 0.015; T2_S = 0.09; k_FS = 8; Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;
figure(1)
subplot(2,3,1)
errorbar(linspace(0,0.25,8),mean(Solution_Bouhrara_MWF(:,:,5),2),std(Solution_Bouhrara_MWF(:,:,5),0,2),'Marker','.','MarkerSize',20,'LineStyle','-','LineWidth',2,'Color',cm(5,:));
hold on;
plot(linspace(0,0.25,8),linspace(0,0.25,8),'LineStyle','--','Color','k','LineWidth',2)
grid on; grid minor; xlim([0 0.25])
get(gca, 'XTick'); set(gca, 'FontSize', 20); get(gca, 'YTick'); set(gca, 'FontSize', 20);
xlabel('MWF','FontSize',20); ylabel('Estimated MWF','FontSize',20)

subplot(2,3,4)
errorbar(linspace(0,0.25,8),(mean(Solution_Bouhrara_MWF(:,:,1),2)/T1_S),(std(Solution_Bouhrara_MWF(:,:,1),0,2)/T1_S),'Marker','.','MarkerSize',20,'LineStyle','-','LineWidth',2,'Color',cm(1,:));
hold on;
errorbar(linspace(0,0.25,8),(mean(Solution_Bouhrara_MWF(:,:,2),2)/T1_F),(std(Solution_Bouhrara_MWF(:,:,2),0,2)/T1_F),'Marker','.','MarkerSize',20,'LineStyle','-','LineWidth',2,'Color',cm(2,:));
errorbar(linspace(0,0.25,8),(mean(Solution_Bouhrara_MWF(:,:,3),2)/T2_S),(std(Solution_Bouhrara_MWF(:,:,3),0,2)/T2_S),'Marker','.','MarkerSize',20,'LineStyle','-','LineWidth',2,'Color',cm(3,:));
errorbar(linspace(0,0.25,8),(mean(Solution_Bouhrara_MWF(:,:,4),2)/T2_F),(std(Solution_Bouhrara_MWF(:,:,4),0,2)/T2_F),'Marker','.','MarkerSize',20,'LineStyle','-','LineWidth',2,'Color',cm(4,:));
errorbar(linspace(0,0.25,8),(mean(Solution_Bouhrara_MWF(:,:,6),2)/k_FS),(std(Solution_Bouhrara_MWF(:,:,6),0,2)/k_FS),'Marker','.','MarkerSize',20,'LineStyle','-','LineWidth',2,'Color',cm(6,:));
plot(linspace(0,0.25,8),linspace(1,1,8),'LineStyle','--','Color','k','LineWidth',2)
grid on; grid minor; xlim([0 0.25])
get(gca, 'XTick'); set(gca, 'FontSize', 20); get(gca, 'YTick'); set(gca, 'FontSize', 20);
xlabel('MWF','FontSize',20); ylabel('Normalised Estimates','FontSize',20)

T1_F = 0.45; T2_F = 0.015; T2_S = 0.09; M0_F = 0.15; k_FS = 8; 
figure(1)
subplot(2,3,2)
errorbar(linspace(1.0,2.0,8),mean(Solution_Bouhrara_T1(:,:,1),2),std(Solution_Bouhrara_T1(:,:,1),0,2),'Marker','.','MarkerSize',20,'LineStyle','-','LineWidth',2,'Color',cm(1,:));
ll = legend('T_{1S}'); ll.AutoUpdate = 'off'; ll.FontSize = 22; legend boxoff; ll.Orientation = 'horizontal'; ll.Position = [0.028282318749481,0.558311319453932,0.051562499254942,0.076023389721474];
hold on;
plot(linspace(1.0,2.0,8),linspace(1.0,2.0,8),'LineStyle','--','Color','k','LineWidth',2)
grid on; grid minor;
get(gca, 'XTick'); set(gca, 'FontSize', 20); get(gca, 'YTick'); set(gca, 'FontSize', 20);
xlabel('T_{1S} (s)','FontSize',20); ylabel('Estimated T_{1S}','FontSize',20)

subplot(2,3,5)
errorbar(linspace(1.0,2.0,8),(mean(Solution_Bouhrara_T1(:,:,2),2)/T1_F),(std(Solution_Bouhrara_T1(:,:,2),0,2)/T1_F),'Marker','.','MarkerSize',20,'LineStyle','-','LineWidth',2,'Color',cm(2,:));
hold on;
errorbar(linspace(1.0,2.0,8),(mean(Solution_Bouhrara_T1(:,:,3),2)/T2_S),(std(Solution_Bouhrara_T1(:,:,3),0,2)/T2_S),'Marker','.','MarkerSize',20,'LineStyle','-','LineWidth',2,'Color',cm(3,:));
errorbar(linspace(1.0,2.0,8),(mean(Solution_Bouhrara_T1(:,:,4),2)/T2_F),(std(Solution_Bouhrara_T1(:,:,4),0,2)/T2_F),'Marker','.','MarkerSize',20,'LineStyle','-','LineWidth',2,'Color',cm(4,:));
errorbar(linspace(1.0,2.0,8),(mean(Solution_Bouhrara_T1(:,:,5),2)/M0_F),(std(Solution_Bouhrara_T1(:,:,5),0,2)/M0_F),'Marker','.','MarkerSize',20,'LineStyle','-','LineWidth',2,'Color',cm(5,:));
errorbar(linspace(1.0,2.0,8),(mean(Solution_Bouhrara_T1(:,:,6),2)/k_FS),(std(Solution_Bouhrara_T1(:,:,6),0,2)/k_FS),'Marker','.','MarkerSize',20,'LineStyle','-','LineWidth',2,'Color',cm(6,:));
plot(linspace(1.0,2.0,8),linspace(1,1,8),'LineStyle','--','Color','k','LineWidth',2)
grid on; grid minor;
get(gca, 'XTick'); set(gca, 'FontSize', 20); get(gca, 'YTick'); set(gca, 'FontSize', 20);
xlabel('T_{1S} (s)','FontSize',20); ylabel('Normalised Estimates','FontSize',20)
ll = legend('T_{1F}','T_{2S}','T_{2F}','MWF','k_{FS}'); ll.FontSize = 22; legend boxoff; ll.Orientation = 'vertical'; ll.Position = [-0.0524,0.2709,0.2266,0.3406]; %ylim([0.8 2])

T1_F = 0.45; T1_S = 1.4; T2_F = 0.015; M0_F = 0.15; k_FS = 8;  
figure(1)
subplot(2,3,3)
errorbar(linspace(0.04,0.1,8),mean(Solution_Bouhrara_T2(:,:,3),2),std(Solution_Bouhrara_T2(:,:,3),0,2),'Marker','.','MarkerSize',20,'LineStyle','-','LineWidth',2,'Color',cm(3,:));
hold on;
plot(linspace(0.04,0.1,8),linspace(0.04,0.1,8),'LineStyle','--','Color','k','LineWidth',2)
grid on; grid minor;
get(gca, 'XTick'); set(gca, 'FontSize', 20); get(gca, 'YTick'); set(gca, 'FontSize', 20);
xlabel('T_{2S} (s)','FontSize',20); ylabel('Estimated T_{2S}','FontSize',20)

subplot(2,3,6)
errorbar(linspace(0.04,0.1,8),(mean(Solution_Bouhrara_T2(:,:,1),2)/T1_S),(std(Solution_Bouhrara_T2(:,:,1),0,2)/T1_S),'Marker','.','MarkerSize',20,'LineStyle','-','LineWidth',2,'Color',cm(1,:));
hold on;
errorbar(linspace(0.04,0.1,8),(mean(Solution_Bouhrara_T2(:,:,2),2)/T1_F),(std(Solution_Bouhrara_T2(:,:,2),0,2)/T1_F),'Marker','.','MarkerSize',20,'LineStyle','-','LineWidth',2,'Color',cm(2,:));
errorbar(linspace(0.04,0.1,8),(mean(Solution_Bouhrara_T2(:,:,4),2)/T2_F),(std(Solution_Bouhrara_T2(:,:,4),0,2)/T2_F),'Marker','.','MarkerSize',20,'LineStyle','-','LineWidth',2,'Color',cm(4,:));
errorbar(linspace(0.04,0.1,8),(mean(Solution_Bouhrara_T2(:,:,5),2)/M0_F),(std(Solution_Bouhrara_T2(:,:,5),0,2)/M0_F),'Marker','.','MarkerSize',20,'LineStyle','-','LineWidth',2,'Color',cm(5,:));
errorbar(linspace(0.04,0.1,8),(mean(Solution_Bouhrara_T2(:,:,6),2)/k_FS),(std(Solution_Bouhrara_T2(:,:,6),0,2)/k_FS),'Marker','.','MarkerSize',20,'LineStyle','-','LineWidth',2,'Color',cm(6,:));
plot(linspace(0.04,0.1,8),linspace(1,1,8),'LineStyle','--','Color','k','LineWidth',2)
grid on; grid minor;
get(gca, 'XTick'); set(gca, 'FontSize', 20); get(gca, 'YTick'); set(gca, 'FontSize', 20);
xlabel('T_{2S} (s)','FontSize',20); ylabel('Normalised Estimates','FontSize',20)
