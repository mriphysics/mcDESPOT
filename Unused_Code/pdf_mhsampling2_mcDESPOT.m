%% Performs MC simulations and generates heat maps to depict cost-function space. (MaxAngle is in degrees and phase-cycling is in radians).

close all; clear all;

% Tissue and sequence parameters.
T1_S = 0.965; %T1S_Fix = [0.7, 0.965, 1.3, 1.7, 2, 2.5]; %T1S_Fix = [0.9, 0.965, 1.1, 1.3, 1.4, 1.5];
T1_F = 0.465; %T1F_Fix = [0.2, 0.25, 0.3, 0.4, 0.465, 0.5]; %T1F_Fix = [0.3, 0.4, 0.465, 0.6, 0.7, 0.8]; 
T2_S = 0.090; %T2S_Fix = [0.075, 0.09, 0.125, 0.15, 0.175, 0.2]; %T2S_Fix = [0.04, 0.06, 0.09, 0.11, 0.13, 0.15]; 
T2_F = 0.012; %T2F_Fix = [0.002, 0.012, 0.02, 0.03, 0.04, 0.045]; %T2F_Fix = [0.01, 0.012, 0.015, 0.02, 0.025, 0.03];
k_FS = 8.000; %kFS_Fix = [0.5, 5, 8, 10, 15, 20]; %kFS_Fix = [(1/0.6), 5, 8, 10, 20, 40];
M0_F = 0.200; 
Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;

% Others are all very similar - chosen to be in middle of Toby's bounds.
% T1_S = 1.2; T1_F = 0.55; T2_S = 0.095; T2_F = 0.02; k_FS = 10; M0_F = 0.175; Delta = pi/8; PC1 = 0 + Delta; PC2 = pi + Delta;

% Rui's acquisition parameters.
TR_SPGR = 7e-3; TR_SSFP = 7e-3; FA_SPGR = deg2rad([6 8 10 12 14 16]); FA_SSFP180 = deg2rad([15 25 35 45 55 65]); FA_SSFP0 = deg2rad([25 55]);

% Deoni's acquisition parameters (2015, SRC Paper)
% TR_SPGR = 5.6e-3; TR_SSFP = 4.4e-3; FA_SPGR = deg2rad([4 5 6 7 9 11 14 18]); FA_SSFP0 = deg2rad([12 16 19 23 27 34 50 70]); FA_SSFP180 = deg2rad([12 16 19 23 27 34 50 70]); 

% Ground-truth signals for mcDESPOT.
SPGR_Data = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'M0_F',M0_F,'k_FS',k_FS);
SSFP_Data_0 = SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);
SSFP_Data_180 = SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);

SSFP_Data = [SSFP_Data_0 ; SSFP_Data_180];

SNR = 30; Sigma = mean(SPGR_Data)/SNR;
SPGR_Data_Noisy = zeros(length(SPGR_Data),1);
for mm = 1:length(SPGR_Data)
    SPGR_Data_Noisy(mm) = SPGR_Data(mm) + (normrnd(0,Sigma));
end
SSFP_Data_Noisy = zeros(length(SSFP_Data),1);
for nn = 1:length(SSFP_Data)
    SSFP_Data_Noisy(nn) = SSFP_Data(nn) + (normrnd(0,Sigma));
end

Data_N = [SPGR_Data_Noisy ; SSFP_Data_Noisy];


%% Perform stochastic region contraction.

Realisations = 100; Trials = 5000; Iterations = 7; N = 50; Runs = 1; Params = 7;
SNR = 30; Sigma = mean(SPGR_Data)/SNR; % Depends on param set.

Solution = zeros(1, Realisations,Params);

delete(gcp('nocreate')); c = parcluster('local'); c.NumWorkers = 8; parpool(c, c.NumWorkers);

tic
for ii = 1:1
    
    for tt = 1:Realisations
        disp(['Realisation: ', num2str(tt)])
        
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
        Data_NN = [SPGR_Data_Norm ; SSFP_Data_Norm];
        
        %SPGR_Data_Norm = SPGR_Data./mean(SPGR_Data);
        %SSFP_Data_Norm = SSFP_Data./mean(SSFP_Data);
        %Data_Noiseless = [SPGR_Data_Norm ; SSFP_Data_Norm];
        
        % Post-normalisation.
        [Solution(ii,tt,:), ~, ~] = SRC_Sim(Trials, Iterations, N, Runs, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_NN);
        
    end

end
toc

%% Direct sampling of cost-function in 2D.

SNR = 30; Sigma = mean(SPGR_Data)/SNR; % Depends on param set.

Upper = [1.5 0.8 0.35 40 0.15 0.03 pi]; Lower = [0.9 0.3 0.001 (1/0.6) 0.04 0.01 -pi]; % Wood.
% Upper = [2.5 0.5 0.30 20 0.2 0.045 pi]; Lower = [0.7 0.2 1e-7 0.5 0.075 0.002 -pi]; % Zhang.
% Upper = [3.5 0.65 0.45 40 0.2 0.06 pi]; Lower = [0.01 0.01 0 1 0.06 0.001 -pi]; % Bouhrara.
% Upper = [5 0.65 0.35 40 0.165 0.03 2*pi]; Lower = [0.9 0.3 0 (1/0.6) 0.05 0.001 0]; % Deoni.

% Pre-normalisation.
Data_O = [SPGR_Data ; SSFP_Data_0 ; SSFP_Data_180];

Steps = 500; nTrials = 1;

T1S_Vector = linspace(Lower(1),Upper(1),Steps); T1F_Vector = linspace(Lower(2),Upper(2),Steps);
M0F_Vector = linspace(Lower(3),Upper(3),Steps); kFS_Vector = linspace(Lower(4),Upper(4),Steps);
T2S_Vector = linspace(Lower(5),Upper(5),Steps); T2F_Vector = linspace(Lower(6),Upper(6),Steps);
Delta_Vector = linspace(Lower(7),Upper(7),Steps);

P_T1 = zeros(Steps,Steps,nTrials); P_T2 = zeros(Steps,Steps,nTrials);
P_T1F = zeros(Steps,Steps,nTrials); P_T1S = zeros(Steps,Steps,nTrials);
P_T2F = zeros(Steps,Steps,nTrials); P_T2S = zeros(Steps,Steps,nTrials);
P_kFS = zeros(Steps,Steps,nTrials);
%P_All = zeros(nTrials,1);

T1S_Rand = zeros(nTrials,1); T1F_Rand = zeros(nTrials,1);
T2S_Rand = zeros(nTrials,1); T2F_Rand = zeros(nTrials,1);
M0F_Rand = zeros(nTrials,1); kFS_Rand = zeros(nTrials,1);
Delta_Rand = zeros(nTrials,1);

for ii = 1:Steps
    disp(['Outer Step Number: ', num2str(ii), '/', num2str(Steps), '.'])
    for jj = 1:Steps
        disp(['     Inner Step Number: ', num2str(jj), '/', num2str(Steps), '.'])
        tic
        for nn = 1:nTrials
            
            %display(['nTrial: ', num2str(nn)])
            % Choose GT as first combination.
            if nn == 1
                T1F_Rand = T1_F; T1S_Rand = T1_S; T2F_Rand = T2_F; T2S_Rand = T2_S; kFS_Rand = k_FS; M0F_Rand = M0_F; Delta_Rand = Delta;
                %T1F_Rand = Solution(2); T1S_Rand = Solution(1); T2F_Rand = Solution(4); T2S_Rand = Solution(3); kFS_Rand = Solution(6); M0F_Rand = Solution(5); Delta_Rand = Solution(7);
                else
                
                % Draw random parameter values from a uniform distribution.
                T1S_Rand(nn) = (Upper(1) - Lower(1)) .* rand(1,1) + Lower(1);
                T1F_Rand(nn) = (Upper(2) - Lower(2)) .* rand(1,1) + Lower(2);
                M0F_Rand(nn) = (Upper(3) - Lower(3)) .* rand(1,1) + Lower(3);
                kFS_Rand(nn) = (Upper(4) - Lower(4)) .* rand(1,1) + Lower(4);
                T2S_Rand(nn) = (Upper(5) - Lower(5)) .* rand(1,1) + Lower(5);
                T2F_Rand(nn) = (Upper(6) - Lower(6)) .* rand(1,1) + Lower(6);
                Delta_Rand(nn) = (Upper(7) - Lower(7)) .* rand(1,1) + Lower(7);
                
            end
            
            P_T1(ii,jj,nn) = (logpdf_mcDESPOT([T1S_Vector(ii) T1F_Vector(jj) M0F_Rand kFS_Rand T2S_Rand T2F_Rand Delta_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
            P_T2(ii,jj,nn) = (logpdf_mcDESPOT([T1S_Rand T1F_Rand M0F_Rand kFS_Rand T2S_Vector(ii) T2F_Vector(jj) Delta_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
            P_T1F(ii,jj,nn) = (logpdf_mcDESPOT([T1S_Rand T1F_Vector(ii) M0F_Vector(jj) kFS_Rand T2S_Rand T2F_Rand Delta_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
            P_T1S(ii,jj,nn) = (logpdf_mcDESPOT([T1S_Vector(ii) T1F_Rand M0F_Vector(jj) kFS_Rand T2S_Rand T2F_Rand Delta_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
            P_T2F(ii,jj,nn) = (logpdf_mcDESPOT([T1S_Rand T1F_Rand M0F_Vector(ii) kFS_Rand T2S_Rand T2F_Vector(jj) Delta_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
            P_T2S(ii,jj,nn) = (logpdf_mcDESPOT([T1S_Rand T1F_Rand M0F_Vector(ii) kFS_Rand T2S_Vector(jj) T2F_Rand Delta_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
            P_kFS(ii,jj,nn) = (logpdf_mcDESPOT([T1S_Rand T1F_Rand M0F_Vector(ii) kFS_Vector(jj) T2S_Rand T2F_Rand Delta_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
            
        end
    end
    
    %P_All(nn) = (logpdf_mcDESPOT([T1S_Rand(nn) T1F_Rand(nn) M0F_Rand(nn) kFS_Rand(nn) T2S_Rand(nn) T2F_Rand(nn) Delta_Rand(nn)], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_N, Sigma));
    
end

%Exp_P_All = exp(P_All);

%% Signal curve analysis.
% PlottingNo = 100;
% [Values, Indexes] = sort(Exp_P_All,'ascend');
% Idx_Picked = Indexes(nTrials-(PlottingNo-1):nTrials);
% T1F_Picked = T1F_Rand(Idx_Picked); T1S_Picked = T1S_Rand(Idx_Picked);
% T2F_Picked = T2F_Rand(Idx_Picked); T2S_Picked = T2S_Rand(Idx_Picked);
% kFS_Picked = kFS_Rand(Idx_Picked); M0F_Picked = M0F_Rand(Idx_Picked);
% Delta_Picked = Delta_Rand(Idx_Picked);
% 
% Idx_Picked2 = Indexes(nTrials-9:nTrials);
% T1F_Picked2 = T1F_Rand(Idx_Picked2); T1S_Picked2 = T1S_Rand(Idx_Picked2);
% T2F_Picked2 = T2F_Rand(Idx_Picked2); T2S_Picked2 = T2S_Rand(Idx_Picked2);
% kFS_Picked2 = kFS_Rand(Idx_Picked2); M0F_Picked2 = M0F_Rand(Idx_Picked2);
% Delta_Picked2 = Delta_Rand(Idx_Picked2);
% 
% figure(1);
% plot(FA_SPGR,SPGR_Data_Noisy/mean(SPGR_Data_Noisy), 'bo','Linewidth',2); hold on
% plot(FA_SSFP0, SSFP0_Data_Noisy/mean(SSFP0_Data_Noisy), 'ko','LineWidth',2)
% plot(FA_SSFP180, SSFP180_Data_Noisy/mean(SSFP180_Data_Noisy), 'ro','LineWidth',2)
% xlabel('FA [rad]', 'FontSize', 12); ylabel('Normalised Signal (a.u.)', 'FontSize',12); 
% ll = legend('SPGR','bSSFP_{0}','bSSFP_{180}'); ll.FontSize = 14; ll.AutoUpdate = 'off'; legend('boxoff');
% 
% SPGR_FM = zeros(length(FA_SPGR),length(Idx_Picked));
% SSFP0_FM = zeros(length(FA_SSFP0),length(Idx_Picked));
% SSFP180_FM = zeros(length(FA_SSFP180),length(Idx_Picked));
% 
% for ii = 1:length(Idx_Picked)
%     
%     SPGR_FM(:,ii) = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1S_Picked(ii),'T1_F',T1F_Picked(ii),'M0_F',M0F_Picked(ii),'k_FS',kFS_Picked(ii));
%     SSFP0_FM(:,ii) = SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1,'T1_S',T1S_Picked(ii),'T2_S',T2S_Picked(ii),'T1_F',T1F_Picked(ii),'T2_F',T2F_Picked(ii),'M0_F',M0F_Picked(ii),'k_FS',kFS_Picked(ii));
%     SSFP180_FM(:,ii) = SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2,'T1_S',T1S_Picked(ii),'T2_S',T2S_Picked(ii),'T1_F',T1F_Picked(ii),'T2_F',T2F_Picked(ii),'M0_F',M0F_Picked(ii),'k_FS',kFS_Picked(ii));
%     
% end
% 
% LineColours = viridis(PlottingNo); colormap viridis; CBar = colorbar; caxis([Values(nTrials-(PlottingNo-1)) Values(nTrials)]); 
% ylabel(CBar,'Likelihood','FontSize',12)
% 
% for jj = 1:length(Idx_Picked)
%    
%     figure(1); plot(FA_SPGR, SPGR_FM(:,jj)/mean(SPGR_FM(:,jj)), '--', 'Color', LineColours(jj,:), 'LineWidth', 1.5);
%     plot(FA_SSFP0, SSFP0_FM(:,jj)/mean(SSFP0_FM(:,jj)), '--', 'Color', LineColours(jj,:), 'LineWidth',0.1);
%     plot(FA_SSFP180, SSFP180_FM(:,jj)/mean(SSFP180_FM(:,jj)), '--', 'Color', LineColours(jj,:), 'LineWidth',0.1); hold on; 
% 
% end
% 
% figure(2);
% subplot(3,2,1); histogram(T1F_Picked,15,'FaceColor','k','FaceAlpha',0.5); xlabel('T_{1F} (s)','FontSize',12); ylabel('Count','FontSize',12); hold on
%                 histogram(T1F_Picked2,5,'FaceColor','r','FaceAlpha',0.5); xlabel('T_{1F} (s)','FontSize',12); ylabel('Count','FontSize',12); hold on
% subplot(3,2,2); histogram(T1S_Picked,15,'FaceColor','k','FaceAlpha',0.5); xlabel('T_{1S} (s)','FontSize',12); ylabel('Count','FontSize',12); hold on
%                 histogram(T1S_Picked2,5,'FaceColor','r','FaceAlpha',0.5); xlabel('T_{1S} (s)','FontSize',12); ylabel('Count','FontSize',12); hold on
% subplot(3,2,3); histogram(T2F_Picked,15,'FaceColor','k','FaceAlpha',0.5); xlabel('T_{2F} (s)','FontSize',12); ylabel('Count','FontSize',12); hold on
%                 histogram(T2F_Picked2,5,'FaceColor','r','FaceAlpha',0.5); xlabel('T_{2F} (s)','FontSize',12); ylabel('Count','FontSize',12); hold on
% subplot(3,2,4); histogram(T2S_Picked,15,'FaceColor','k','FaceAlpha',0.5); xlabel('T_{2S} (s)','FontSize',12); ylabel('Count','FontSize',12); hold on
%                 histogram(T2S_Picked2,5,'FaceColor','r','FaceAlpha',0.5); xlabel('T_{2S} (s)','FontSize',12); ylabel('Count','FontSize',12); hold on
% subplot(3,2,5); histogram(M0F_Picked,15,'FaceColor','k','FaceAlpha',0.5); xlabel('M_{0F}','FontSize',12); ylabel('Count','FontSize',12); hold on
%                 histogram(M0F_Picked2,5,'FaceColor','r','FaceAlpha',0.5); xlabel('M_{0F}','FontSize',12); ylabel('Count','FontSize',12); hold on
% subplot(3,2,6); histogram(kFS_Picked,15,'FaceColor','k','FaceAlpha',0.5); xlabel('k_{FS} (s^{-1})','FontSize',12); ylabel('Count','FontSize',12); hold on
%                 histogram(kFS_Picked2,3,'FaceColor','r','FaceAlpha',0.5); xlabel('k_{FS} (s^{-1})','FontSize',12); ylabel('Count','FontSize',12); hold on
% ll = legend('Top 100', 'Top 10'); ll.FontSize = 14;

%% Heat map analysis and plots.

LineSettings = ['w', 'w', 'w', 'w', 'w', 'w', 'w'];

% Bounds Order: T1S_LB, T1S_UB, T1F_LB, T1F_UB, T2S_LB, T2S_UB, T2F_LB, T2F_UB, M0F_LB, M0F_UB, kFS_LB, kFS_UB;
%T1S_LB = Bounds(:,1); T1S_UB = Bounds(:,2); T1F_LB = Bounds(:,3); T1F_UB = Bounds(:,4);
%T2S_LB = Bounds(:,5); T2S_UB = Bounds(:,6); T2F_LB = Bounds(:,7); T2F_UB = Bounds(:,8);
%M0F_LB = Bounds(:,9); M0F_UB = Bounds(:,10); kFS_LB = Bounds(:,11); kFS_UB = Bounds(:,12);

figure(1); subplot(3,3,1);
MeanCF_T1 = mean(exp(P_T1),3); MaxCF_T1 = max(exp(P_T1),[],3);
imagesc([min(T1F_Vector) max(T1F_Vector)],[min(T1S_Vector) max(T1S_Vector)],MaxCF_T1); colormap(viridis); colorbar; shading interp; xlabel('T_{1F} (s)','FontSize',12); ylabel('T_{1S} (s)','FontSize',12); hold on;
% for rr = 1:length(T1S_LB)
%     rectangle('Position',[T1F_LB(rr) T1S_LB(rr) (T1F_UB(rr)-T1F_LB(rr)) (T1S_UB(rr)-T1S_LB(rr))],'EdgeColor',LineSettings(rr),'LineWidth',1)
% end
% for ss = 1:size(Subsets,1)
%     plot(Subsets(ss,2),Subsets(ss,1),'m+','LineWidth',1)
% end
%plot(Solution(:,2), Solution(:,1),'m.', 'MarkerSize', 2); 
plot(T1_F,T1_S,'w.', 'MarkerSize', 20);

subplot(3,3,2);
MeanCF_T2 = mean(exp(P_T2),3); MaxCF_T2 = max(exp(P_T2),[],3);
imagesc([min(T2F_Vector) max(T2F_Vector)],[min(T2S_Vector) max(T2S_Vector)],MaxCF_T2); colormap(viridis); colorbar; shading interp; xlabel('T_{2F} (s)','FontSize',12); ylabel('T_{2S} (s)','FontSize',12); hold on;
% for rr = 1:length(T1S_LB)
%     rectangle('Position',[T2F_LB(rr) T2S_LB(rr) (T2F_UB(rr)-T2F_LB(rr)) (T2S_UB(rr)-T2S_LB(rr))],'EdgeColor',LineSettings(rr),'LineWidth',1)
% end
% for ss = 1:size(Subsets,1)
%     plot(Subsets(ss,4),Subsets(ss,3),'m+','LineWidth',1)
% end
%plot(Solution(:,4), Solution(:,3),'m.', 'MarkerSize', 2); 
plot(T2_F,T2_S,'w.', 'MarkerSize', 20);

subplot(3,3,3)
MeanCF_kFS = mean(exp(P_kFS),3); MaxCF_kFS = max(exp(P_kFS),[],3);
imagesc([min(kFS_Vector) max(kFS_Vector)],[min(M0F_Vector) max(M0F_Vector)],MaxCF_kFS); colormap(viridis); colorbar; shading interp; xlabel('k_{FS} (s^{-1})','FontSize',12); ylabel('M_{0F}','FontSize',12); hold on;
% for rr = 1:length(T1S_LB)
%     rectangle('Position',[kFS_LB(rr) M0F_LB(rr) (kFS_UB(rr)-kFS_LB(rr)) (M0F_UB(rr)-M0F_LB(rr))],'EdgeColor',LineSettings(rr),'LineWidth',1)
% end
% for ss = 1:size(Subsets,1)
%     plot(Subsets(ss,6),Subsets(ss,5),'m+','LineWidth',1)
% end
%plot(Solution(:,6), Solution(:,5),'m.', 'MarkerSize', 2); 
plot(k_FS,M0_F,'w.', 'MarkerSize', 20);

subplot(3,3,4);
MeanCF_T1F = mean(exp(P_T1F),3); MaxCF_T1F = max(exp(P_T1F),[],3);
imagesc([min(M0F_Vector) max(M0F_Vector)],[min(T1F_Vector) max(T1F_Vector)],MaxCF_T1F); colormap(viridis); colorbar; shading interp; xlabel('M_{0F}','FontSize',12); ylabel('T_{1F} (s)','FontSize',12); hold on;
% for rr = 1:length(T1S_LB)
%     rectangle('Position',[M0F_LB(rr) T1F_LB(rr) (M0F_UB(rr)-M0F_LB(rr)) (T1F_UB(rr)-T1F_LB(rr))],'EdgeColor',LineSettings(rr),'LineWidth',1)
% end
% for ss = 1:size(Subsets,1)
%     plot(Subsets(ss,5),Subsets(ss,2),'m+','LineWidth',1)
% end
%plot(Solution(:,5), Solution(:,2),'m.', 'MarkerSize', 2); 
plot(M0_F,T1_F,'w.', 'MarkerSize', 20);

subplot(3,3,5);
MeanCF_T1S = mean(exp(P_T1S),3); MaxCF_T1S = max(exp(P_T1S),[],3);
imagesc([min(M0F_Vector) max(M0F_Vector)],[min(T1S_Vector) max(T1S_Vector)],MaxCF_T1S); colormap(viridis); colorbar; shading interp; xlabel('M_{0F}','FontSize',12); ylabel('T_{1S} (s)','FontSize',12); hold on;
% for rr = 1:length(T1S_LB)
%     rectangle('Position',[M0F_LB(rr) T1S_LB(rr) (M0F_UB(rr)-M0F_LB(rr)) (T1S_UB(rr)-T1S_LB(rr))],'EdgeColor',LineSettings(rr),'LineWidth',1)
% end
% for ss = 1:size(Subsets,1)
%     plot(Subsets(ss,5),Subsets(ss,1),'m+','LineWidth',1)
% end
%plot(Solution(:,5), Solution(:,1),'m.', 'MarkerSize', 2); 
plot(M0_F,T1_S,'w.', 'MarkerSize', 20);

subplot(3,3,6);
MeanCF_T2F = mean(exp(P_T2F),3); MaxCF_T2F = max(exp(P_T2F),[],3);
imagesc([min(T2F_Vector) max(T2F_Vector)],[min(M0F_Vector) max(M0F_Vector)],MaxCF_T2F); colormap(viridis); colorbar; shading interp; xlabel('T_{2F} (s)','FontSize',12); ylabel('M_{0F}','FontSize',12); hold on;
% for rr = 1:length(T1S_LB)
%     rectangle('Position',[T2F_LB(rr) M0F_LB(rr) (T2F_UB(rr)-T2F_LB(rr)) (M0F_UB(rr)-M0F_LB(rr))],'EdgeColor',LineSettings(rr),'LineWidth',1)
% end
% for ss = 1:size(Subsets,1)
%     plot(Subsets(ss,4),Subsets(ss,5),'m+','LineWidth',1)
% end
%plot(Solution(:,4), Solution(:,5),'m.', 'MarkerSize', 2); 
plot(T2_F,M0_F,'w.', 'MarkerSize', 20);

subplot(3,3,8);
MeanCF_T2S = mean(exp(P_T2S),3); MaxCF_T2S = max(exp(P_T2S),[],3);
imagesc([min(T2S_Vector) max(T2S_Vector)],[min(M0F_Vector) max(M0F_Vector)],MaxCF_T2S); colormap(viridis); colorbar; shading interp; xlabel('T_{2S} (s)','FontSize',12); ylabel('M_{0F}','FontSize',12); hold on;
% for rr = 1:length(T1S_LB)
%     rectangle('Position',[T2S_LB(rr) M0F_LB(rr) (T2S_UB(rr)-T2S_LB(rr)) (M0F_UB(rr)-M0F_LB(rr))],'EdgeColor',LineSettings(rr),'LineWidth',1)
% end
% for ss = 1:size(Subsets,1)
%     plot(Subsets(ss,3),Subsets(ss,5),'m+','LineWidth',1)
% end
%plot(Solution(:,3), Solution(:,5),'m.', 'MarkerSize', 2); 
plot(T2_S,M0_F,'w.', 'MarkerSize', 20);

%% Bound progression plots.

% figure(Realisations+tt)
% subplot(4,1,1)
% plot(1:Iterations+1, T1F_LB, 'k--','LineWidth',2)
% hold on;
% plot(1:Iterations+1, T1F_UB, 'k--','LineWidth',2)
% plot(1:Iterations+1, T1S_LB, 'r--','LineWidth',2)
% plot(1:Iterations+1, T1S_UB, 'r--','LineWidth',2)
% xlabel('Iterations', 'FontSize', 12); ylabel('T_{1}', 'FontSize', 12)
% ll = legend('T_{1F}','','T_{1S}',''); ll.FontSize = 12;
% hline(T1_F,'k'); hline(T1_S,'r')
%
% subplot(4,1,2)
% plot(1:Iterations+1, T2F_LB, 'k--','LineWidth',2)
% hold on;
% plot(1:Iterations+1, T2F_UB, 'k--','LineWidth',2)
% plot(1:Iterations+1, T2S_LB, 'r--','LineWidth',2)
% plot(1:Iterations+1, T2S_UB, 'r--','LineWidth',2)
% xlabel('Iterations', 'FontSize', 12); ylabel('T_{2}', 'FontSize', 12)
% ll = legend('T_{2F}','','T_{2S}',''); ll.FontSize = 12;
% hline(T2_F,'k'); hline(T2_S,'r')
%
% subplot(4,1,3)
% plot(1:Iterations+1, M0F_LB, 'k--','LineWidth',2)
% hold on;
% plot(1:Iterations+1, M0F_UB, 'k--','LineWidth',2)
% xlabel('Iterations', 'FontSize', 12); ylabel('M_{0F}', 'FontSize', 12)
% hline(M0_F,'k')
%
% subplot(4,1,4)
% plot(1:Iterations+1, kFS_LB, 'k--','LineWidth',2)
% hold on;
% plot(1:Iterations+1, kFS_UB, 'k--','LineWidth',2)
% xlabel('Iterations', 'FontSize', 12); ylabel('k_{FS}', 'FontSize', 12)
% hline(k_FS,'k')


%% Parameter estimate histograms.

% Solution_ReOrder = Solution;
% 
% Solution_ReOrder(3,:,:) = Solution(1,:,:); Solution_ReOrder(1,:,:) = Solution(3,:,:);
% 
% for pp = 1:3
%     
%     figure(2); FaceColours = ['k','m','c'];
%     subplot(2,3,1); histogram(Solution_ReOrder(pp,:,1), 'BinWidth', (T1_S/60), 'EdgeColor' ,'k', 'FaceColor', FaceColours(pp), 'FaceAlpha', 0.6); xlabel('T_{1S} (s)','FontSize',12); ylabel('Count','FontSize',12); hold on; ll = legend('200','50','30'); ll.FontSize = 12; vline(T1_S, 'r--');
%     subplot(2,3,2); histogram(Solution_ReOrder(pp,:,2), 'BinWidth', (T1_F/60), 'EdgeColor' ,'k', 'FaceColor', FaceColours(pp), 'FaceAlpha', 0.6); xlabel('T_{1F} (s)','FontSize',12); ylabel('Count','FontSize',12); hold on; ll = legend('200','50','30'); ll.FontSize = 12; vline(T1_F, 'r--');
%     subplot(2,3,3); histogram(Solution_ReOrder(pp,:,3), 'BinWidth', (T2_S/60), 'EdgeColor' ,'k', 'FaceColor', FaceColours(pp), 'FaceAlpha', 0.6); xlabel('T_{2S} (s)','FontSize',12); ylabel('Count','FontSize',12); hold on; ll = legend('200','50','30'); ll.FontSize = 12; vline(T2_S, 'r--');
%     subplot(2,3,4); histogram(Solution_ReOrder(pp,:,4), 'BinWidth', (T2_F/30), 'EdgeColor' ,'k', 'FaceColor', FaceColours(pp), 'FaceAlpha', 0.6); xlabel('T_{2F} (s)','FontSize',12); ylabel('Count','FontSize',12); hold on; ll = legend('200','50','30'); ll.FontSize = 12; vline(T2_F, 'r--');
%     subplot(2,3,5); histogram(Solution_ReOrder(pp,:,5), 'BinWidth', (M0_F/50), 'EdgeColor' ,'k', 'FaceColor', FaceColours(pp), 'FaceAlpha', 0.6); xlabel('M_{0F}','FontSize',12); ylabel('Count','FontSize',12); hold on; ll = legend('200','50','30'); ll.FontSize = 12; vline(M0_F, 'r--');
%     subplot(2,3,6); histogram(Solution_ReOrder(pp,:,6), 'BinWidth', (k_FS/20), 'EdgeColor' ,'k', 'FaceColor', FaceColours(pp), 'FaceAlpha', 0.6); xlabel('k_{FS} (s^{-1})','FontSize',12); ylabel('Count','FontSize',12); hold on; ll = legend('200','50','30'); ll.FontSize = 12; vline(k_FS, 'r--');
%     
% end

%% MWF box plots.
%BP_Vector = [M0F_All ; M0F_FixedT1F ; M0F_FixedT1S ; M0F_FixedT2F ; M0F_FixedT2S ; M0F_FixedkFS ; M0F_FixedWT1F ; M0F_FixedWT1S ; M0F_FixedWT2F ; M0F_FixedWT2S ; M0F_FixedWkFS ; M0F_FixedT1FT2F ; M0F_FixedT1ST2S ; M0F_FixedT1FT1S ; M0F_FixedT2FT2S];
%Plotting_Vector = [ones(100,1); 2*ones(100,1); 3*ones(100,1); 4*ones(100,1); 5*ones(100,1); 6*ones(100,1); 7*ones(100,1); 8*ones(100,1); 9*ones(100,1); 10*ones(100,1); 11*ones(100,1); 12*ones(100,1); 13*ones(100,1); 14*ones(100,1); 15*ones(100,1)];
%boxplot(BP_Vector,Plotting_Vector)
%Labels = {'All','T1F','T1S','T2F','T2S','kFS','2T1F','2T1S','2T2F','2T2S','2kFS','T1F&T2F','T1S&T2S','T1F&T1S','T2F&T2S'};
%set(gca,'XTickLabel',Labels)
%hold on; hline(0.2, 'm--');

% For Wood's bounds.
% figure(1);
% subplot(2,3,1); boxplot([Solution(:,5), SolutionT1F(1,:,5)', SolutionT1F(2,:,5)', SolutionT1F(3,:,5)', SolutionT1F(4,:,5)', SolutionT1F(5,:,5)', SolutionT1F(6,:,5)']);
% a = get(get(gca,'children'),'children'); t = get(a,'tag'); box = a(18); set(box,'Color','k')
% ylabel('MWF'); Labels = {'All','0.3','0.4','0.465','0.6','0.7','0.8'}; set(gca,'XTickLabel',Labels); xlabel('T_{1F} (s)'); hline(0.2, 'm--');
% subplot(2,3,2); boxplot([Solution(:,5), SolutionT1S(1,:,5)', SolutionT1S(2,:,5)', SolutionT1S(3,:,5)', SolutionT1S(4,:,5)', SolutionT1S(5,:,5)', SolutionT1S(6,:,5)']); 
% a = get(get(gca,'children'),'children'); box = a(19); set(box,'Color','k')
% ylabel('MWF'); Labels = {'All','0.9','0.965','1.1','1.3','1.4','1.5'}; set(gca,'XTickLabel',Labels); xlabel('T_{1S} (s)'); hline(0.2, 'm--');
% subplot(2,3,3); boxplot([Solution(:,5), SolutionT2F(1,:,5)', SolutionT2F(2,:,5)', SolutionT2F(3,:,5)', SolutionT2F(4,:,5)', SolutionT2F(5,:,5)', SolutionT2F(6,:,5)']); 
% a = get(get(gca,'children'),'children'); box = a(19); set(box,'Color','k')
% ylabel('MWF'); Labels = {'All','10','12','15','20','25','30'}; set(gca,'XTickLabel',Labels); xlabel('T_{2F} (ms)'); hline(0.2, 'm--');
% subplot(2,3,4); boxplot([Solution(:,5), SolutionT2S(1,:,5)', SolutionT2S(2,:,5)', SolutionT2S(3,:,5)', SolutionT2S(4,:,5)', SolutionT2S(5,:,5)', SolutionT2S(6,:,5)']); 
% a = get(get(gca,'children'),'children'); box = a(18); set(box,'Color','k')
% ylabel('MWF'); Labels = {'All','40','60','90','110','130','150'}; set(gca,'XTickLabel',Labels); xlabel('T_{2S} (ms)'); hline(0.2, 'm--');
% subplot(2,3,5); boxplot([Solution(:,5), SolutionkFS(1,:,5)', SolutionkFS(2,:,5)', SolutionkFS(3,:,5)', SolutionkFS(4,:,5)', SolutionkFS(5,:,5)', SolutionkFS(6,:,5)']); 
% a = get(get(gca,'children'),'children'); box = a(18); set(box,'Color','k')
% ylabel('MWF'); Labels = {'All','(1/0.6)','5','8','10','20','40'}; set(gca,'XTickLabel',Labels); xlabel('k_{FS} (s^{-1})'); hline(0.2, 'm--');
% 
% % For Zhang's bounds.
% figure(2);
% subplot(2,3,1); boxplot([Solution(:,5), SolutionT1F(1,:,5)', SolutionT1F(2,:,5)', SolutionT1F(3,:,5)', SolutionT1F(4,:,5)', SolutionT1F(5,:,5)', SolutionT1F(6,:,5)']);
% a = get(get(gca,'children'),'children'); box = a(16); set(box,'Color','k')
% ylabel('MWF'); Labels = {'All','0.2','0.25','0.3','0.4','0.465','0.5'}; set(gca,'XTickLabel',Labels); xlabel('T_{1F} (s)'); hline(0.2, 'm--');
% subplot(2,3,2); boxplot([Solution(:,5), SolutionT1S(1,:,5)', SolutionT1S(2,:,5)', SolutionT1S(3,:,5)', SolutionT1S(4,:,5)', SolutionT1S(5,:,5)', SolutionT1S(6,:,5)']); 
% a = get(get(gca,'children'),'children'); box = a(19); set(box,'Color','k')
% ylabel('MWF'); Labels = {'All','0.7','0.965','1.3','1.7','2.0','2.5'}; set(gca,'XTickLabel',Labels); xlabel('T_{1S} (s)'); hline(0.2, 'm--');
% subplot(2,3,3); boxplot([Solution(:,5), SolutionT2F(1,:,5)', SolutionT2F(2,:,5)', SolutionT2F(3,:,5)', SolutionT2F(4,:,5)', SolutionT2F(5,:,5)', SolutionT2F(6,:,5)']); 
% a = get(get(gca,'children'),'children'); box = a(19); set(box,'Color','k')
% ylabel('MWF'); Labels = {'All','2','12','20','30','40','45'}; set(gca,'XTickLabel',Labels); xlabel('T_{2F} (ms)'); hline(0.2, 'm--');
% subplot(2,3,4); boxplot([Solution(:,5), SolutionT2S(1,:,5)', SolutionT2S(2,:,5)', SolutionT2S(3,:,5)', SolutionT2S(4,:,5)', SolutionT2S(5,:,5)', SolutionT2S(6,:,5)']); 
% a = get(get(gca,'children'),'children'); box = a(19); set(box,'Color','k')
% ylabel('MWF'); Labels = {'All','75','90','125','150','175','200'}; set(gca,'XTickLabel',Labels); xlabel('T_{2S} (ms)'); hline(0.2, 'm--');
% subplot(2,3,5); boxplot([Solution(:,5), SolutionkFS(1,:,5)', SolutionkFS(2,:,5)', SolutionkFS(3,:,5)', SolutionkFS(4,:,5)', SolutionkFS(5,:,5)', SolutionkFS(6,:,5)']); 
% a = get(get(gca,'children'),'children'); box = a(18); set(box,'Color','k')
% ylabel('MWF'); Labels = {'All','0.5','5','8','10','15','20'}; set(gca,'XTickLabel',Labels); xlabel('k_{FS} (s^{-1})'); hline(0.2, 'm--');


%%%%% OUTTAKES %%%%%

%TR_SPGR = 5e-3; TR_SSFP = 5e-3;
%FA_SPGR = deg2rad([3 4 5 6 7 9 13 18]); FA_SSFP0 = deg2rad([10 13 17 20 23 30 43 60]); FA_SSFP180 = deg2rad([10 13 17 20 23 30 43 60]); % Assumed: M0_S = 0.8; k_SF = (M0_F*k_FS)/M0_S;
% TR_SPGR = 5.2e-3; TR_SSFP = 5.2e-3;
% FA_SPGR = [3, 4, 5, 6, 7, 9, 13, 18]; FA_SSFP0 = [2, 5, 10, 15, 20, 30, 40, 50]; FA_SSFP180 = [2, 5, 10, 15, 20, 30, 40, 50];
% T1_S = 1.15; T1_F = 0.4; T2_S = 0.08; T2_F = 0.02; k_FS = 9; M0_F = [0.1, 0.2, 0.3];
% Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;
%FA_SSFP0 = deg2rad([12 16 21 27 33 41 52.5 65]); FA_SSFP180 = deg2rad([12 16 21 27 33 41 52.5 65]); %2PCP.
%PC3 = pi/2 + Delta; PC4 = 3*pi/2 + Delta; FA_SSFP0 = deg2rad([12 21 33 52.5]); FA_SSFP180 = deg2rad([12 21 33 52.5]); FA_SSFP90 = deg2rad([16 27 41 65]); FA_SSFP270 = deg2rad([16 27 41 65]); %4PCP.
%FA_SSFP180 = deg2rad([12 12 16 16 21 21 27 27 33 33 41 41 52.5 52.5 65 65]); %1PCP.

%SSFP_Data_90 = SSFP_SteadyState(FA_SSFP90, TR_SSFP, PC3,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);
%SSFP_Data_270 = SSFP_SteadyState(FA_SSFP270, TR_SSFP, PC4,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);
%SSFP_Data_90_Norm = SSFP_Data_90./mean(SSFP_Data_90); SSFP_Data_270_Norm = SSFP_Data_270./mean(SSFP_Data_270);
% SSFP_Data_90_Norm ; SSFP_Data_270_Norm
% FA_SSFP90, FA_SSFP270,