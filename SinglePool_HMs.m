%%% Generates single-pool heat maps for Supplementary Figure 4.

close all; clear all;

% Tissue and sequence parameters.
T1 = 1; T2 = 0.1; M0 = 1; Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;
TR_SPGR = 6.5e-3; TR_SSFP = 6.5e-3; FA_SPGR = deg2rad([2 4 6 8 10 12 14 16 18 20]); FA_SSFP180 = deg2rad([2 6 14 22 30 38 46 54 62 70]); FA_SSFP0 = deg2rad([2 6 14 22 30 38 46 54 62 70]);

% Ground-truth signals for mcDESPOT.
SPGR_Data = SPGR_SP_SteadyState(FA_SPGR, TR_SPGR,'T1',T1,'M0',M0);
SSFP_Data_0 = SSFP_SP_SteadyState(FA_SSFP0, TR_SSFP, PC1,'T1',T1,'T2',T2,'M0',M0);
SSFP_Data_180 = SSFP_SP_SteadyState(FA_SSFP180, TR_SSFP, PC2,'T1',T1,'T2',T2,'M0',M0);

% Concatenate SSFP signals.
SSFP_Data = [SSFP_Data_0 ; SSFP_Data_180];

%% Direct sampling of cost-function in 2D. Modify for signal fit plots.
SNR = 30; Sigma = mean(SPGR_Data)/SNR;

Upper = [1.1 1.1 0.11]; Lower = [0.9 0.9 0.09];

% Pre-normalisation.
Data_O = [SPGR_Data ; SSFP_Data_0 ; SSFP_Data_180];

Steps = 49; nTrials = 100000;

T1_Vector = linspace(Lower(1),Upper(1),Steps); M0_Vector = linspace(Lower(2),Upper(2),Steps); 
T2_Vector = linspace(Lower(3),Upper(3),Steps);

P_T1 = zeros(Steps,Steps,nTrials); P_T2 = zeros(Steps,Steps,nTrials); 
P_T1T2 = zeros(Steps,Steps,nTrials);

for ii = 1:Steps
    disp(['Outer Step Number: ', num2str(ii), '/', num2str(Steps), '.'])
    for jj = 1:Steps
        parfor nn = 1:nTrials

            % Draw random parameter values from a uniform distribution.
            T1_Rand = (Upper(1) - Lower(1)) .* rand(1,1) + Lower(1);
            M0_Rand = (Upper(2) - Lower(2)) .* rand(1,1) + Lower(2);
            T2_Rand = (Upper(3) - Lower(3)) .* rand(1,1) + Lower(3);
                            
            P_T1(ii,jj,nn) = (logpdf_SinglePool([T1_Vector(ii) M0_Vector(jj) T2_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
            P_T2(ii,jj,nn) = (logpdf_SinglePool([T1_Rand M0_Vector(ii) T2_Vector(jj)], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
            P_T1T2(ii,jj,nn) = (logpdf_SinglePool([T1_Vector(ii) M0_Rand T2_Vector(jj)], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
            
        end
    end
end

%% Heat map analysis and plots.

coloraxis = ([0 1]);

figure(3); subplot(1,3,1); %MaxCF_T1 = max(exp(P_T1),[],3);
imagesc([min(M0_Vector) max(M0_Vector)],[min(T1_Vector) max(T1_Vector)],MaxCF_T1_NoiseInd,coloraxis); shading interp; xlabel('M_{0}','FontSize',16); ylabel('T_{1} (s)','FontSize',16); hold on;
plot(M0,T1,'.','Color',[0 0.5 0], 'MarkerSize', 20); get(gca, 'XTick'); set(gca, 'FontSize', 18); get(gca, 'YTick'); set(gca, 'FontSize', 18); xlim([0.95 1.05])
axis square

subplot(1,3,2); %MaxCF_T2 = max(exp(P_T2),[],3);
imagesc([min(T2_Vector) max(T2_Vector)],[min(M0_Vector) max(M0_Vector)],MaxCF_T2_NoiseInd,coloraxis); shading interp; xlabel('T_{2} (s)','FontSize',16); ylabel('M_{0}','FontSize',16); hold on;
plot(T2,M0,'.','Color',[0 0.5 0], 'MarkerSize', 20); get(gca, 'XTick'); set(gca, 'FontSize', 18); get(gca, 'YTick'); set(gca, 'FontSize', 18); xlim([0.095 0.105])
axis square

subplot(1,3,3); %MaxCF_T1T2 = max(exp(P_T1T2),[],3);
imagesc([min(T2_Vector) max(T2_Vector)],[min(T1_Vector) max(T1_Vector)],MaxCF_T1T2_NoiseInd,coloraxis); shading interp; xlabel('T_{2} (s)','FontSize',16); ylabel('T_{1} (s)','FontSize',16); hold on;
plot(T2,T1,'.','Color',[0 0.5 0], 'MarkerSize', 20); get(gca, 'XTick'); set(gca, 'FontSize', 18); get(gca, 'YTick'); set(gca, 'FontSize', 18); xlim([0.095 0.105])
axis square

colormap(flipud(magma)); 
hcb = colorbar('Position',[0.909583090080568,0.291713961407491,0.013774574153008,0.450624290578888],'FontSize',14); hcb.FontSize = 18;
colorTitleHandle = get(hcb,'Title'); titleString = '\eta';
set(colorTitleHandle,'String',titleString,'FontSize',30); 
