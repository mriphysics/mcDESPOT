clear all
close all

TR_SPGR = 5.2e-3; TR_SSFP = 5.2e-3; % Make different?
FA_SPGR = deg2rad([3 4 5 6 7 9 13 18]); 
FA_SSFP = deg2rad([2 5 10 15 20 30 40 50]);
T1_S = 1.15; T1_F = 0.4; T1_B = 1;
T2_S = 0.08; T2_F = 0.02;
M0_F = 0.25; M0_S = 0.55; M0_B = 0.2;
k_FS = 9; k_SB = 5; k_FB = 5; k_SF = (M0_F*k_FS)/M0_S;

R1_F = 1/T1_F; R1_B = 1/T1_B; R1_S = 1/T1_S; dMzB = 0;
G = 1.4e-5; Gamma = 2 * pi * 42.57747892e6; B1 = 13e-6;
T_RF = deg2rad(50)/(Gamma * B1); W = (pi/T_RF) * (Gamma * B1)^2 * T_RF * G * (T_RF/TR_SPGR);

kSF_Prediction = (M0_F*k_FS)/M0_S + (M0_F*k_FB*k_SB)/(M0_B*W + M0_F*k_FB + M0_S*k_SB + M0_B*R1_B);
kFS_Prediction = k_FS + (M0_S*k_FB*k_SB)/(M0_B*W + M0_F*k_FB + M0_S*k_SB + M0_B*R1_B);
R1F_Prediction = (M0_B*R1_B*R1_F + M0_B*R1_F*W + M0_B*R1_B*k_FB + M0_F*R1_F*k_FB + M0_S*R1_F*k_SB + M0_B*W*k_FB)/(M0_B*W + M0_F*k_FB + M0_S*k_SB + M0_B*R1_B);
R1S_Prediction = (M0_B*R1_B*R1_S + M0_B*R1_S*W + M0_F*R1_S*k_FB + M0_B*R1_B*k_SB + M0_S*R1_S*k_SB + M0_B*W*k_SB)/(M0_B*W + M0_F*k_FB + M0_S*k_SB + M0_B*R1_B);
M0F_Prediction = (M0_F*(M0_B*R1_B*R1_F - dMzB*k_FB + M0_B*R1_F*W + M0_B*R1_B*k_FB + M0_F*R1_F*k_FB + M0_S*R1_F*k_SB))/(M0_B*R1_B*R1_F + M0_B*R1_F*W + M0_B*R1_B*k_FB + M0_F*R1_F*k_FB + M0_S*R1_F*k_SB + M0_B*W*k_FB);
M0S_Prediction = (M0_S*(M0_B*R1_B*R1_S - dMzB*k_SB + M0_B*R1_S*W + M0_F*R1_S*k_FB + M0_B*R1_B*k_SB + M0_S*R1_S*k_SB))/(M0_B*R1_B*R1_S + M0_B*R1_S*W + M0_F*R1_S*k_FB + M0_B*R1_B*k_SB + M0_S*R1_S*k_SB + M0_B*W*k_SB);
T1F_Prediction = 1/R1F_Prediction;
T1S_Prediction = 1/R1S_Prediction;

% GT-signals.
SPGR_Data = TRF_SPGR_steady_state_MT_M0(FA_SPGR,TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP_Data = TRF_SSFP_steady_state_MT_M0(FA_SSFP,TR_SSFP,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP_Data_180 = TRF_SSFP_steady_state_180_MT_M0(FA_SSFP,TR_SSFP,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
Data = [SPGR_Data ; SSFP_Data; SSFP_Data_180];

Sigma = 1/300; %SNR = 500.

%Upper = [1.5 0.75 0.5 0.75 20 10 0.2 0.1];
%Lower = [0 0 0 0 0 0 0 0];
Upper = [1.5 0.8 0.35 0.75 40 40 0.15 0.03];
Lower = [0.9 0.3 0.001 0.001 (1/0.6) (1/0.6) 0.04 0.01];

%% Direct sampling of cost-function in 2D.

Steps = 100;

T1S_Vector = linspace(Lower(1),Upper(1),Steps);
T1F_Vector = linspace(Lower(2),Upper(2),Steps);
M0F_Vector = linspace(Lower(3),Upper(3),Steps);
kFS_Vector = linspace(Lower(5),Upper(5),Steps);
kSF_Vector = linspace(Lower(6),Upper(6),Steps);
T2S_Vector = linspace(Lower(7),Upper(7),Steps);
T2F_Vector = linspace(Lower(8),Upper(8),Steps);

P_T1S = zeros(length(T1S_Vector),length(T1S_Vector));
P_T1F = zeros(length(T1S_Vector),length(T1S_Vector));
P_kFS = zeros(length(T1S_Vector),length(T1S_Vector));
P_kSF = zeros(length(T1S_Vector),length(T1S_Vector));
P_T2S = zeros(length(T1S_Vector),length(T1S_Vector));
P_T2F = zeros(length(T1S_Vector),length(T1S_Vector));

disp (['Estimated Run Time: ' num2str(0.0109*Steps^2) ' seconds']);

for ii = 1:Steps
    tic
    for jj = 1:Steps
        
        P_T1S(ii,jj) = (logpdf([T1S_Vector(ii) T1_F M0F_Vector(jj) M0_S k_FS k_SF T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        P_T1F(ii,jj) = (logpdf([T1_S T1F_Vector(ii) M0F_Vector(jj) M0_S k_FS k_SF T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        P_kFS(ii,jj) = (logpdf([T1_S T1_F M0F_Vector(ii) M0_S kFS_Vector(jj) k_SF T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        P_kSF(ii,jj) = (logpdf([T1_S T1_F M0F_Vector(ii) M0_S k_FS kSF_Vector(jj) T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        P_T2S(ii,jj) = (logpdf([T1_S T1_F M0F_Vector(ii) M0_S k_FS k_SF T2S_Vector(jj) T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        P_T2F(ii,jj) = (logpdf([T1_S T1_F M0F_Vector(ii) M0_S k_FS k_SF T2_S T2F_Vector(jj)], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        
    end
   toc 
end

Matrix_kFS = exp(P_kFS); Matrix_kSF = exp(P_kSF);
Matrix_T1F = exp(P_T1F); Matrix_T1S = exp(P_T1S);
Matrix_T2F = exp(P_T2F); Matrix_T2S = exp(P_T2S);
Threshold = 0.999999;

[RowX_kFS, ColX_kFS] = find(Matrix_kFS > Threshold);
[RowX_kSF, ColX_kSF] = find(Matrix_kSF > Threshold);
[RowX_T1F, ColX_T1F] = find(Matrix_T1F > Threshold);
[RowX_T1S, ColX_T1S] = find(Matrix_T1S > Threshold);
[RowX_T2F, ColX_T2F] = find(Matrix_T2F > Threshold);
[RowX_T2S, ColX_T2S] = find(Matrix_T2S > Threshold);

figure(5)
subplot(3,2,1)
imagesc([min(M0F_Vector) max(M0F_Vector)],[min(T1S_Vector) max(T1S_Vector)],exp(P_T1S)); colorbar; shading interp; xlabel('M_{0F}'); ylabel('T_{1S}'); hold on; plot(M0_F,T1_S,'w.', 'MarkerSize', 20); plot(M0F_Prediction,T1S_Prediction,'k.', 'MarkerSize', 20);
hold on; scatter(M0F_Vector(ColX_T1S), T1S_Vector(RowX_T1S), 'rx')

subplot(3,2,2)
imagesc([min(M0F_Vector) max(M0F_Vector)],[min(T1F_Vector) max(T1F_Vector)],exp(P_T1F)); colorbar; shading interp; xlabel('M_{0F}'); ylabel('T_{1F}'); hold on; plot(M0_F,T1_F,'w.', 'MarkerSize', 20); plot(M0F_Prediction,T1F_Prediction,'k.', 'MarkerSize', 20);
hold on; scatter(M0F_Vector(ColX_T1F), T1F_Vector(RowX_T1F), 'rx')

subplot(3,2,3)
imagesc([min(kFS_Vector) max(kFS_Vector)],[min(M0F_Vector) max(M0F_Vector)],exp(P_kFS)); colorbar; shading interp; xlabel('k_{FS}'); ylabel('M_{0F}'); hold on; plot(k_FS,M0_F,'w.', 'MarkerSize', 20); plot(kFS_Prediction,M0F_Prediction,'k.', 'MarkerSize', 20);
hold on; scatter(kFS_Vector(ColX_kFS), M0F_Vector(RowX_kFS), 'rx')

subplot(3,2,4)
imagesc([min(kSF_Vector) max(kSF_Vector)],[min(M0F_Vector) max(M0F_Vector)],exp(P_kSF)); colorbar; shading interp; xlabel('k_{SF}'); ylabel('M_{0F}'); hold on; plot(k_SF,M0_F,'w.', 'MarkerSize', 20); plot(kSF_Prediction,M0F_Prediction,'k.', 'MarkerSize', 20);
hold on; scatter(kSF_Vector(ColX_kSF), M0F_Vector(RowX_kSF), 'rx')

subplot(3,2,5)
imagesc([min(T2S_Vector) max(T2S_Vector)],[min(M0F_Vector) max(M0F_Vector)],exp(P_T2S)); colorbar; shading interp; xlabel('T_{2S}'); ylabel('M_{0F}'); hold on; plot(T2_S,M0_F,'w.', 'MarkerSize', 20); plot(T2_S,M0F_Prediction,'k.', 'MarkerSize', 20);
hold on; scatter(T2S_Vector(ColX_T2S), M0F_Vector(RowX_T2S), 'rx')

subplot(3,2,6)
imagesc([min(T2F_Vector) max(T2F_Vector)],[min(M0F_Vector) max(M0F_Vector)],exp(P_T2F)); colorbar; shading interp; xlabel('T_{2F}'); ylabel('M_{0F}'); hold on; plot(T2_F,M0_F,'w.', 'MarkerSize', 20); plot(T2_F,M0F_Prediction,'k.', 'MarkerSize', 20);
hold on; scatter(T2F_Vector(ColX_T2F), M0F_Vector(RowX_T2F), 'rx')