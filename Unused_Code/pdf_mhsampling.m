%% Generates heat maps to depict cost-function space. MaxAngle is in degrees and phase-cycling is in radians.

close all; clear all

% Tissue and sequence parameters.
TR_SPGR = 5e-3; TR_SSFP = 5e-3;
FA_SPGR = deg2rad([3 4 5 6 7 9 13 18]); FA_SSFP0 = deg2rad([10 13 17 20 23 30 43 60]); FA_SSFP180 = deg2rad([10 13 17 20 23 30 43 60]); MaxAngle = 60;
T1_S = 0.965; T1_F = 0.465; T1_B = 1; T2_S = 0.09; T2_F = 0.012; M0_F = 0.2; k_FS = 8; M0_S = 0.7; M0_B = 0.1; k_SB = 5; k_FB = 5; k_SF = (M0_F*k_FS)/M0_S; 
PC1 = 0; PC2 = pi;

%% Ground-truth signals. 

% CSMT GT.
SPGR_Data = CS_SPGR_SteadyState(FA_SPGR, TR_SPGR, MaxAngle,'T1_S',T1_S,'T1_F',T1_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP_Data_0 = CS_SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1, MaxAngle,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP_Data_180 = CS_SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2, MaxAngle,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
Data = [SPGR_Data ; SSFP_Data_0 ; SSFP_Data_180];

Sigma = 1/300; %SNR = 30.
Upper = [1.25 0.5 0.30 0.80 20 10 0.15 0.03]; Lower = [0.20 0.1 0.05 0.05 1 1 0.04 0.01];

% Calulcate parameter predictions from derived expressions - only valid for CSMT.
R1_F = 1/T1_F; R1_B = 1/T1_B; R1_S = 1/T1_S; dMzB = 0;
G = 1.4e-5; Gamma = 2 * pi * 42.57747892e6; B1 = 13e-6;
T_RF = deg2rad(MaxAngle)/(Gamma * B1); W = (pi/T_RF) * (Gamma * B1)^2 * T_RF * G * (T_RF/TR_SSFP);
kSF_Prediction = (M0_F*k_FS)/M0_S + (M0_F*k_FB*k_SB)/(M0_B*W + M0_F*k_FB + M0_S*k_SB + M0_B*R1_B);
kFS_Prediction = k_FS + (M0_S*k_FB*k_SB)/(M0_B*W + M0_F*k_FB + M0_S*k_SB + M0_B*R1_B);
R1F_Prediction = (M0_B*R1_B*R1_F + M0_B*R1_F*W + M0_B*R1_B*k_FB + M0_F*R1_F*k_FB + M0_S*R1_F*k_SB + M0_B*W*k_FB)/(M0_B*W + M0_F*k_FB + M0_S*k_SB + M0_B*R1_B);
R1S_Prediction = (M0_B*R1_B*R1_S + M0_B*R1_S*W + M0_F*R1_S*k_FB + M0_B*R1_B*k_SB + M0_S*R1_S*k_SB + M0_B*W*k_SB)/(M0_B*W + M0_F*k_FB + M0_S*k_SB + M0_B*R1_B);
M0F_Prediction = (M0_F*(M0_B*R1_B*R1_F - dMzB*k_FB + M0_B*R1_F*W + M0_B*R1_B*k_FB + M0_F*R1_F*k_FB + M0_S*R1_F*k_SB))/(M0_B*R1_B*R1_F + M0_B*R1_F*W + M0_B*R1_B*k_FB + M0_F*R1_F*k_FB + M0_S*R1_F*k_SB + M0_B*W*k_FB);
M0S_Prediction = (M0_S*(M0_B*R1_B*R1_S - dMzB*k_SB + M0_B*R1_S*W + M0_F*R1_S*k_FB + M0_B*R1_B*k_SB + M0_S*R1_S*k_SB))/(M0_B*R1_B*R1_S + M0_B*R1_S*W + M0_F*R1_S*k_FB + M0_B*R1_B*k_SB + M0_S*R1_S*k_SB + M0_B*W*k_SB);
T1F_Prediction = 1/R1F_Prediction; T1S_Prediction = 1/R1S_Prediction;

%% Direct sampling of cost-function in 2D.

Steps = 100; nTrials = 100;
T1S_Vector = linspace(Lower(1),Upper(1),Steps); T1F_Vector = linspace(Lower(2),Upper(2),Steps);
M0F_Vector = linspace(Lower(3),Upper(3),Steps); M0S_Vector = linspace(Lower(4),Upper(4),Steps);
kFS_Vector = linspace(Lower(5),Upper(5),Steps); kSF_Vector = linspace(Lower(6),Upper(6),Steps);
T2S_Vector = linspace(Lower(7),Upper(7),Steps); T2F_Vector = linspace(Lower(8),Upper(8),Steps);

T1S_Rand = (Upper(1) - Lower(1)) .* rand(nTrials,1) + Lower(1);
T1F_Rand = (Upper(2) - Lower(2)) .* rand(nTrials,1) + Lower(2);
M0F_Rand = (Upper(3) - Lower(3)) .* rand(nTrials,1) + Lower(3);
M0S_Rand = (Upper(4) - Lower(4)) .* rand(nTrials,1) + Lower(4);
kFS_Rand = (Upper(5) - Lower(5)) .* rand(nTrials,1) + Lower(5);
kSF_Rand = (Upper(6) - Lower(6)) .* rand(nTrials,1) + Lower(6);
T2S_Rand = (Upper(7) - Lower(7)) .* rand(nTrials,1) + Lower(7);
T2F_Rand = (Upper(8) - Lower(8)) .* rand(nTrials,1) + Lower(8);

% P_T1S = zeros(length(T1S_Vector),length(T1S_Vector));
% P_T1F = zeros(length(T1S_Vector),length(T1S_Vector));
% P_kFS = zeros(length(T1S_Vector),length(T1S_Vector));
% P_kSF = zeros(length(T1S_Vector),length(T1S_Vector));
% P_T2S = zeros(length(T1S_Vector),length(T1S_Vector));
% P_T2F = zeros(length(T1S_Vector),length(T1S_Vector));
P_T1 = zeros(length(T1S_Vector),length(T1F_Vector),nTrials);
P_T2 = zeros(length(T2S_Vector),length(T2F_Vector),nTrials);
P_EX = zeros(length(kFS_Vector),length(kSF_Vector),nTrials);
P_M0 = zeros(length(M0F_Vector),length(M0S_Vector),nTrials);

for nn = 1:nTrials
    disp(['Trial Number ', num2str(nn)])
    for ii = 1:length(T1S_Vector)
        disp(['Step Number: ', num2str(ii)])
        for jj = 1:length(T1F_Vector)
            tic
            % Use predictions to define slice for CSMT.
            %P_T1S(ii,jj) = (logpdf([T1S_Vector(ii) T1F_Prediction M0F_Vector(jj) M0S_Prediction kFS_Prediction kSF_Prediction T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
            %P_T1F(ii,jj) = (logpdf([T1S_Prediction T1F_Vector(ii) M0F_Vector(jj) M0S_Prediction kFS_Prediction kSF_Prediction T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
            %P_kFS(ii,jj) = (logpdf([T1S_Prediction T1F_Prediction M0F_Vector(ii) M0S_Prediction kFS_Vector(jj) kSF_Prediction T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
            %P_kSF(ii,jj) = (logpdf([T1S_Prediction T1F_Prediction M0F_Vector(ii) M0S_Prediction kFS_Prediction kSF_Vector(jj) T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
            %P_T2S(ii,jj) = (logpdf([T1S_Prediction T1F_Prediction M0F_Vector(ii) M0S_Prediction kFS_Prediction kSF_Prediction T2S_Vector(jj) T2_F], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
            %P_T2F(ii,jj) = (logpdf([T1S_Prediction T1F_Prediction M0F_Vector(ii) M0S_Prediction kFS_Prediction kSF_Prediction T2_S T2F_Vector(jj)], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
            P_T1(ii,jj,nn) = (logpdf([T1S_Vector(ii) T1F_Vector(jj) M0F_Rand(nn) M0S_Rand(nn) kFS_Rand(nn) kSF_Rand(nn) T2S_Rand(nn) T2F_Rand(nn)], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
            P_T2(ii,jj,nn) = (logpdf([T1S_Rand(nn) T1F_Rand(nn) M0F_Rand(nn) M0S_Rand(nn) kFS_Rand(nn) kSF_Rand(nn) T2S_Vector(ii) T2F_Vector(jj)], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
            P_EX(ii,jj,nn) = (logpdf([T1S_Rand(nn) T1F_Rand(nn) M0F_Rand(nn) M0S_Rand(nn) kFS_Vector(ii) kSF_Vector(jj) T2S_Rand(nn) T2F_Rand(nn)], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
            P_M0(ii,jj,nn) = (logpdf([T1S_Rand(nn) T1F_Rand(nn) M0F_Vector(ii) M0S_Vector(jj) kFS_Rand(nn) kSF_Rand(nn) T2S_Rand(nn) T2F_Rand(nn)], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
            toc
        end
    end
end

% Matrix_kFS = exp(P_kFS); Matrix_kSF = exp(P_kSF);
% Matrix_T1F = exp(P_T1F); Matrix_T1S = exp(P_T1S);
% Matrix_T2F = exp(P_T2F); Matrix_T2S = exp(P_T2S);
Matrix_T1 = exp(P_T1); Matrix_T2 = exp(P_T2); Matrix_EX = exp(P_EX); Matrix_M0 = exp(P_M0); 

% figure(5)
% subplot(3,2,1)
% imagesc([min(M0F_Vector) max(M0F_Vector)],[min(T1S_Vector) max(T1S_Vector)],exp(P_T1S)); colorbar; shading interp; xlabel('M_{0F}'); ylabel('T_{1S}'); hold on; plot(M0_F,T1_S,'w.', 'MarkerSize', 20); plot(M0F_Prediction,T1S_Prediction,'k.', 'MarkerSize', 20);
% subplot(3,2,2)
% imagesc([min(M0F_Vector) max(M0F_Vector)],[min(T1F_Vector) max(T1F_Vector)],exp(P_T1F)); colorbar; shading interp; xlabel('M_{0F}'); ylabel('T_{1F}'); hold on; plot(M0_F,T1_F,'w.', 'MarkerSize', 20); plot(M0F_Prediction,T1F_Prediction,'k.', 'MarkerSize', 20);
% subplot(3,2,3)
% imagesc([min(kFS_Vector) max(kFS_Vector)],[min(M0F_Vector) max(M0F_Vector)],exp(P_kFS)); colorbar; shading interp; xlabel('k_{FS}'); ylabel('M_{0F}'); hold on; plot(k_FS,M0_F,'w.', 'MarkerSize', 20); plot(kFS_Prediction,M0F_Prediction,'k.', 'MarkerSize', 20);
% subplot(3,2,4)
% imagesc([min(kSF_Vector) max(kSF_Vector)],[min(M0F_Vector) max(M0F_Vector)],exp(P_kSF)); colorbar; shading interp; xlabel('k_{SF}'); ylabel('M_{0F}'); hold on; plot(k_SF,M0_F,'w.', 'MarkerSize', 20); plot(kSF_Prediction,M0F_Prediction,'k.', 'MarkerSize', 20);
% subplot(3,2,5)
% imagesc([min(T2S_Vector) max(T2S_Vector)],[min(M0F_Vector) max(M0F_Vector)],exp(P_T2S)); colorbar; shading interp; xlabel('T_{2S}'); ylabel('M_{0F}'); hold on; plot(T2_S,M0_F,'w.', 'MarkerSize', 20); plot(T2_S,M0F_Prediction,'k.', 'MarkerSize', 20);
% subplot(3,2,6)
% imagesc([min(T2F_Vector) max(T2F_Vector)],[min(M0F_Vector) max(M0F_Vector)],exp(P_T2F)); colorbar; shading interp; xlabel('T_{2F}'); ylabel('M_{0F}'); hold on; plot(T2_F,M0_F,'w.', 'MarkerSize', 20); plot(T2_F,M0F_Prediction,'k.', 'MarkerSize', 20);

% figure(6)
% subplot(2,2,1)
% imagesc([min(T1F_Vector) max(T1F_Vector)],[min(T1S_Vector) max(T1S_Vector)],exp(P_T1)); colorbar; shading interp; xlabel('T_{1F}'); ylabel('T_{1S}'); hold on; plot(T1_F,T1_S,'w.', 'MarkerSize', 20); plot(T1F_Prediction,T1S_Prediction,'k.', 'MarkerSize', 20);
% subplot(2,2,2)
% imagesc([min(T2F_Vector) max(T2F_Vector)],[min(T2S_Vector) max(T2S_Vector)],exp(P_T2)); colorbar; shading interp; xlabel('T_{2F}'); ylabel('T_{2S}'); hold on; plot(T2_F,T2_S,'w.', 'MarkerSize', 20); plot(T2_F,T2_S,'k.', 'MarkerSize', 20);
% subplot(2,2,3)
% imagesc([min(kSF_Vector) max(kSF_Vector)],[min(kFS_Vector) max(kFS_Vector)],exp(P_EX)); colorbar; shading interp; xlabel('k_{SF}'); ylabel('k_{FS}'); hold on; plot(k_SF,k_FS,'w.', 'MarkerSize', 20); plot(kSF_Prediction,kFS_Prediction,'k.', 'MarkerSize', 20);
% subplot(2,2,4)
% imagesc([min(M0S_Vector) max(M0S_Vector)],[min(M0F_Vector) max(M0F_Vector)],exp(P_M0)); colorbar; shading inte;rp; xlabel('M_{0S}'); ylabel('M_{0F}'); hold on; plot(M0_S,M0_F,'w.', 'MarkerSize', 20); plot(M0S_Prediction,M0F_Prediction,'k.', 'MarkerSize', 20);

%%%% ANALYSIS OF HM STACK GOES HERE %%%
figure(10); subplot(1,2,1);
MeanCF_T1 = mean(exp(P_T1),3); SDCF_T1 = std(exp(P_T1),0,3); MaxCF_T1 = max(exp(P_T1),[],3); 
MaxVal_T1 = max(max(MaxCF_T1)); [Values_T1, Index_T1] = find(MaxCF_T1 > 0.5); Percentage_T1 = (length(Values_T1)/(size(P_T1,1)*size(P_T1,2)))*100;
imagesc([min(T1F_Vector) max(T1F_Vector)],[min(T1S_Vector) max(T1S_Vector)],MaxCF_T1); colormap(magma); colorbar; shading interp; xlabel('T_{1F}','FontSize',12); ylabel('T_{1S}','FontSize',12); hold on; plot(T1_F,T1_S,'w.', 'MarkerSize', 20); plot(T1F_Prediction,T1S_Prediction,'k.', 'MarkerSize', 20);
title('Max Projection of HMs','FontSize',12)
figure(10); subplot(1,2,2);
histogram(MaxCF_T1, 100, 'FaceColor', 'r', 'FaceAlpha',0.5,'EdgeColor','k','EdgeAlpha',0.5);
xlabel('Likelihood','FontSize',12); ylabel('Count','FontSize',12)
title(['Max Value is ',num2str(MaxVal_T1),' and ', num2str(Percentage_T1),'% > 0.5'],'FontSize',12)

figure(11); subplot(1,2,1);
MeanCF_T2 = mean(exp(P_T2),3); SDCF_T2 = std(exp(P_T2),0,3); MaxCF_T2 = max(exp(P_T2),[],3); 
MaxVal_T2 = max(max(MaxCF_T2)); [Values_T2, Index_T2] = find(MaxCF_T2 > 0.5); Percentage_T2 = (length(Values_T2)/(size(P_T2,1)*size(P_T2,2)))*100;
imagesc([min(T2F_Vector) max(T2F_Vector)],[min(T2S_Vector) max(T2S_Vector)],MaxCF_T2); colormap(magma); colorbar; shading interp; xlabel('T_{2F}','FontSize',12); ylabel('T_{2S}','FontSize',12); hold on; plot(T2_F,T2_S,'w.', 'MarkerSize', 20); plot(T2_F,T2_S,'k.', 'MarkerSize', 20);
title('Max Projection of HMs','FontSize',12)
figure(11); subplot(1,2,2);
histogram(MaxCF_T2, 100, 'FaceColor', 'r', 'FaceAlpha',0.5,'EdgeColor','k','EdgeAlpha',0.5);
xlabel('Likelihood','FontSize',12); ylabel('Count','FontSize',12)
title(['Max Value is ',num2str(MaxVal_T2),' and ', num2str(Percentage_T2),'% > 0.5'],'FontSize',12)

figure(12); subplot(1,2,1);
MeanCF_EX = mean(exp(P_EX),3); SDCF_EX = std(exp(P_EX),0,3); MaxCF_EX = max(exp(P_EX),[],3); 
MaxVal_EX = max(max(MaxCF_EX)); [Values_EX, Index_EX] = find(MaxCF_EX > 0.5); Percentage_EX = (length(Values_EX)/(size(P_EX,1)*size(P_EX,2)))*100;
imagesc([min(kSF_Vector) max(kSF_Vector)],[min(kFS_Vector) max(kFS_Vector)],MaxCF_EX); colormap(magma); colorbar; shading interp; xlabel('k_{SF}','FontSize',12); ylabel('k_{FS}','FontSize',12); hold on; plot(k_SF,k_FS,'w.', 'MarkerSize', 20); plot(kSF_Prediction,kFS_Prediction,'k.', 'MarkerSize', 20);
title('Max Projection of HMs','FontSize',12)
figure(12); subplot(1,2,2);
histogram(MaxCF_EX, 100, 'FaceColor', 'r', 'FaceAlpha',0.5,'EdgeColor','k','EdgeAlpha',0.5);
xlabel('Likelihood','FontSize',12); ylabel('Count','FontSize',12)
title(['Max Value is ',num2str(MaxVal_EX),' and ', num2str(Percentage_EX),'% > 0.5'],'FontSize',12)

figure(13); subplot(1,2,1);
MeanCF_M0 = mean(exp(P_M0),3); SDCF_M0 = std(exp(P_M0),0,3); MaxCF_M0 = max(exp(P_M0),[],3); 
MaxVal_M0 = max(max(MaxCF_M0)); [Values_M0, Index_M0] = find(MaxCF_M0 > 0.5); Percentage_M0 = (length(Values_M0)/(size(P_M0,1)*size(P_M0,2)))*100;
imagesc([min(M0S_Vector) max(M0S_Vector)],[min(M0F_Vector) max(M0F_Vector)],MaxCF_M0); colormap(magma); colorbar; shading interp; xlabel('M_{0S}','FontSize',12); ylabel('M_{0F}','FontSize',12); hold on; plot(M0_S,M0_F,'w.', 'MarkerSize', 20); plot(M0S_Prediction,M0F_Prediction,'k.', 'MarkerSize', 20);
title('Max Projection of HMs','FontSize',12)
figure(13); subplot(1,2,2);
histogram(MaxCF_M0, 100, 'FaceColor', 'r', 'FaceAlpha',0.5,'EdgeColor','k','EdgeAlpha',0.5);
xlabel('Likelihood','FontSize',12); ylabel('Count','FontSize',12)
title(['Max Value is ',num2str(MaxVal_M0),' and ', num2str(Percentage_M0),'% > 0.5'],'FontSize',12)

%% Other signal models.

%SPGR_Data = B1_SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
%SSFP_Data_0 = B1_SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
%SSFP_Data_180 = B1_SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
%SSFP_Data_90 = B1_SSFP_SteadyState(FA_SSFP90, TR_SSFP, PC3,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
%SSFP_Data_270 = B1_SSFP_SteadyState(FA_SSFP270, TR_SSFP, PC4,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
%SPGR_Data = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'M0_F',M0_F,'k_FS',k_FS);
%SSFP_Data_0 = SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);
%SSFP_Data_180 = SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);
%SSFP_Data_90 = SSFP_SteadyState(FA_SSFP90, TR_SSFP, PC3,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SF',k_SF);
%SSFP_Data_270 = SSFP_SteadyState(FA_SSFP270, TR_SSFP, PC4,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SF',k_SF);
%SPGR_Data = SPGR_SP_SteadyState(FA_SPGR, TR_SPGR,'T1',T1,'M0',M0);
%SSFP_Data_0 = SSFP_SP_SteadyState(FA_SSFP0, TR_SSFP, PC1,'T1',T1,'T2',T2,'M0',M0);
%SSFP_Data_180 = SSFP_SP_SteadyState(FA_SSFP180, TR_SSFP, PC2,'T1',T1,'T2',T2,'M0',M0);
%SSFP_Data_90 = SSFP_SP_SteadyState(FA_SSFP90, TR_SSFP, PC3,'T1',T1,'T2',T2,'M0',M0);
%SSFP_Data_270 = SSFP_SP_SteadyState(FA_SSFP270, TR_SSFP, PC4,'T1',T1,'T2',T2,'M0',M0);

%% Function call for every other case but CSMT.

%P_T1S(ii,jj) = (logpdf([T1S_Vector(ii) T1_F M0F_Vector(jj) k_FS T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
%P_T1F(ii,jj) = (logpdf([T1_S T1F_Vector(ii) M0F_Vector(jj) k_FS T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
%P_kFS(ii,jj) = (logpdf([T1_S T1_F M0F_Vector(ii) kFS_Vector(jj) T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
%P_kSF(ii,jj) = (logpdf([T1_S T1_F M0F_Vector(ii) k_FS T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
%P_T2S(ii,jj) = (logpdf([T1_S T1_F M0F_Vector(ii) k_FS T2S_Vector(jj) T2_F], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
%P_T2F(ii,jj) = (logpdf([T1_S T1_F M0F_Vector(ii) k_FS T2_S T2F_Vector(jj)], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
%P_T1(ii,jj) = (logpdf([T1S_Vector(ii) T1_F M0_F k_FS T2S_Vector(jj) T2_F], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
%P_T2(ii,jj) = (logpdf([T1_S T1_F M0_F M0_S k_FS k_SF T2S_Vector(ii) T2F_Vector(jj)], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
%P_EX(ii,jj) = (logpdf([T1_S T1_F M0_F M0_S kFS_Vector(ii) kSF_Vector(jj) T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
%P_M0(ii,jj) = (logpdf([T1_S T1_F M0F_Vector(jj) M0S_Vector(ii) k_FS k_SF T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));

%% Thesholding to pick out highest likelihood points.

% Threshold = 0.99999;
% [RowX_kFS, ColX_kFS] = find(Matrix_kFS > Threshold);
% [RowX_kSF, ColX_kSF] = find(Matrix_kSF > Threshold);
% [RowX_T1F, ColX_T1F] = find(Matrix_T1F > Threshold);
% [RowX_T1S, ColX_T1S] = find(Matrix_T1S > Threshold);
% [RowX_T2F, ColX_T2F] = find(Matrix_T2F > Threshold);
% [RowX_T2S, ColX_T2S] = find(Matrix_T2S > Threshold);
%hold on; scatter(M0F_Vector(ColX_T1S), T1S_Vector(RowX_T1S), 'rx')
%hold on; scatter(M0F_Vector(ColX_T1F), T1F_Vector(RowX_T1F), 'rx')
%hold on; scatter(kFS_Vector(ColX_kFS), M0F_Vector(RowX_kFS), 'rx')
%hold on; scatter(kSF_Vector(ColX_kSF), M0F_Vector(RowX_kSF), 'rx')
%hold on; scatter(T2S_Vector(ColX_T2S), M0F_Vector(RowX_T2S), 'rx')
%hold on; scatter(T2F_Vector(ColX_T2F), M0F_Vector(RowX_T2F), 'rx')

%% Plotting for '3D HM'.

% figure(6)
% subplot(3,2,1); 
% imagesc([min(T2S_Vector) max(T2S_Vector)],[min(T1F_Vector) max(T2S_Vector)],exp(P_T1(:,:,1))); colormap(magma); colorbar; shading interp; xlabel('T_{1F}'); ylabel('T_{2S}'); hold on; plot(T1_F,T2_S,'w.', 'MarkerSize', 20); title('M_{0F} = 0.05');
% subplot(3,2,2); 
% imagesc([min(T2S_Vector) max(T2S_Vector)],[min(T1F_Vector) max(T2S_Vector)],exp(P_T1(:,:,2))); colormap(magma); colorbar; shading interp; xlabel('T_{1F}'); ylabel('T_{2S}'); hold on; plot(T1_F,T2_S,'w.', 'MarkerSize', 20); title('M_{0F} = 0.10');
% subplot(3,2,3); 
% imagesc([min(T2S_Vector) max(T2S_Vector)],[min(T1F_Vector) max(T2S_Vector)],exp(P_T1(:,:,3))); colormap(magma); colorbar; shading interp; xlabel('T_{1F}'); ylabel('T_{2S}'); hold on; plot(T1_F,T2_S,'w.', 'MarkerSize', 20); title('M_{0F} = 0.15');
% subplot(3,2,4); 
% imagesc([min(T2S_Vector) max(T2S_Vector)],[min(T1F_Vector) max(T2S_Vector)],exp(P_T1(:,:,4))); colormap(magma); colorbar; shading interp; xlabel('T_{1F}'); ylabel('T_{2S}'); hold on; plot(T1_F,T2_S,'w.', 'MarkerSize', 20); title('M_{0F} = 0.20');
% subplot(3,2,5); 
% imagesc([min(T2S_Vector) max(T2S_Vector)],[min(T2S_Vector) max(T2S_Vector)],exp(P_T1(:,:,5))); colormap(magma); colorbar; shading interp; xlabel('T_{1F}'); ylabel('T_{2S}'); hold on; plot(T1_F,T2_S,'w.', 'MarkerSize', 20); title('M_{0F} = 0.25');

%% Single-pool code.

% T1_Vector = linspace(Lower(1),Upper(1),Steps);
% M0_Vector = linspace(Lower(2),Upper(2),Steps);
% T2_Vector = linspace(Lower(3),Upper(3),Steps);
% P_T1M0 = zeros(length(T1_Vector),length(T1_Vector));
% P_T1T2 = zeros(length(T1_Vector),length(T1_Vector));
% P_M0T2 = zeros(length(T1_Vector),length(T1_Vector));
% 
% % CHANGE INPUTS DEPENDING ON PHASE-CYCLING PATTERN.
% for ii = 1:length(T1_Vector)
%     tic
%     for jj = 1:length(T1_Vector)
%         P_T1M0(ii,jj) = (logpdf([T1_Vector(ii) M0_Vector(jj) T2], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
%         P_T1T2(ii,jj) = (logpdf([T1_Vector(ii) M0 T2_Vector(jj)], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
%         P_M0T2(ii,jj) = (logpdf([T1 M0_Vector(ii) T2_Vector(jj)], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
%     end
%     toc
% end
% 
% Matrix_T1M0 = exp(P_T1M0); Matrix_T1T2 = exp(P_T1T2); Matrix_M0T2 = exp(P_M0T2);
% 
% figure(5)
% subplot(3,1,1)
% imagesc([min(M0_Vector) max(M0_Vector)],[min(T1_Vector) max(T1_Vector)],exp(P_T1M0)); colorbar; shading interp; xlabel('M_{0}'); ylabel('T_{1}');
% subplot(3,1,2)
% imagesc([min(T2_Vector) max(T2_Vector)],[min(T1_Vector) max(T1_Vector)],exp(P_T1T2)); colorbar; shading interp; xlabel('T_{2}'); ylabel('T_{1}');
% subplot(3,1,3)
% imagesc([min(T2_Vector) max(T2_Vector)],[min(M0_Vector) max(M0_Vector)],exp(P_M0T2)); colorbar; shading interp; xlabel('T_{2}'); ylabel('M_{0}');

%% Direct sampling of cost-function in 1D.

% T1S_Vector = linspace(0.2,0.3,1000); T1F_Vector = linspace(0.14,0.22,1000);
% M0F_Vector = linspace(0.08,0.16,1000); M0S_Vector = linspace(0.1,0.18,1000);
% kFS_Vector = linspace(10,11,1000); kSF_Vector = linspace(4.3,5.1,1000);
% T2S_Vector = linspace(0.06,0.1,1000); T2F_Vector = linspace(0.01,0.03,1000);
% 
% P_T1S = zeros(length(T1S_Vector),1); P_T1F = zeros(length(T1S_Vector),1);
% P_M0F = zeros(length(T1S_Vector),1); P_M0S = zeros(length(T1S_Vector),1);
% P_kFS = zeros(length(T1S_Vector),1); P_kSF = zeros(length(T1S_Vector),1);
% P_T2S = zeros(length(T1S_Vector),1); P_T2F = zeros(length(T1S_Vector),1);
% 
% tic
% for ii = 1:length(T1S_Vector)
% 
%     P_T1S(ii) = (logpdf([T1S_Vector(ii) T1F_Prediction M0F_Prediction M0S_Prediction kFS_Prediction kSF_Prediction T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
%     P_T1F(ii) = (logpdf([T1S_Prediction T1F_Vector(ii) M0F_Prediction M0S_Prediction kFS_Prediction kSF_Prediction T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
%     P_M0F(ii) = (logpdf([T1S_Prediction T1F_Prediction M0F_Vector(ii) M0S_Prediction kFS_Prediction kSF_Prediction T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
%     P_M0S(ii) = (logpdf([T1S_Prediction T1F_Prediction M0F_Prediction M0S_Vector(ii) kFS_Prediction kSF_Prediction T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
%     P_kFS(ii) = (logpdf([T1S_Prediction T1F_Prediction M0F_Prediction M0S_Prediction kFS_Vector(ii) kSF_Prediction T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
%     P_kSF(ii) = (logpdf([T1S_Prediction T1F_Prediction M0F_Prediction M0S_Prediction kFS_Prediction kSF_Vector(ii) T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
%     P_T2S(ii) = (logpdf([T1S_Prediction T1F_Prediction M0F_Prediction M0S_Prediction kFS_Prediction kSF_Prediction T2S_Vector(ii) T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
%     P_T2F(ii) = (logpdf([T1S_Prediction T1F_Prediction M0F_Prediction M0S_Prediction kFS_Prediction kSF_Prediction T2_S T2F_Vector(ii)], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
% 
% end
% toc
% 
% figure(4);
% subplot(2,4,1); plot(T1S_Vector,exp(P_T1S),'k'); xlabel('T_{1S}'); ylabel('PDF Value'); hold on; vline(T1S_Prediction, 'r', 'Pred')
% subplot(2,4,2); plot(T1F_Vector,exp(P_T1F),'k'); xlabel('T_{1F}'); ylabel('PDF Value'); hold on; vline(T1F_Prediction, 'r', 'Pred')
% subplot(2,4,4); plot(M0F_Vector,exp(P_M0F),'k'); xlabel('M_{0F}'); ylabel('PDF Value'); hold on; vline(M0F_Prediction, 'r', 'Pred')
% subplot(2,4,3); plot(M0S_Vector,exp(P_M0S),'k'); xlabel('M_{0S}'); ylabel('PDF Value'); hold on; vline(M0S_Prediction, 'r', 'Pred')
% subplot(2,4,5); plot(kFS_Vector,exp(P_kFS),'k'); xlabel('k_{FS}'); ylabel('PDF Value'); hold on; vline(kFS_Prediction, 'r', 'Pred')
% subplot(2,4,6); plot(kSF_Vector,exp(P_kSF),'k'); xlabel('k_{SF}'); ylabel('PDF Value'); hold on; vline(kSF_Prediction, 'r', 'Pred')
% subplot(2,4,7); plot(T2S_Vector,exp(P_T2S),'k'); xlabel('T_{2S}'); ylabel('PDF Value'); hold on; vline(T2_S, 'r', 'Pred')
% subplot(2,4,8); plot(T2F_Vector,exp(P_T2F),'k'); xlabel('T_{2F}'); ylabel('PDF Value'); hold on; vline(T2_F, 'r', 'Pred')
