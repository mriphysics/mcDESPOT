% When running, update signal function to calculate (not read-in) kSF and M0S.

clear all
close all

TR_SPGR = 5.2e-3; TR_SSFP = 5.2e-3;
FA_SPGR = deg2rad([3 4 5 6 7 9 13 18]); 
FA_SSFP = deg2rad([2 5 10 15 20 30 40 50]);

%TR_SPGR = 6.5e-3; TR_SSFP = 5e-3; FA_SPGR = deg2rad([2 4 6 8 10 12 14]); FA_SSFP = deg2rad([6 14 22 30 38 46 54 62 70]);

%T1_S = 1.15; T1_F = 0.4; T2_S = 0.08; T2_F = 0.02;
%M0_F = 0.25; k_FS = 9; 
dMzB = 0;

T1_S = 0.965; T1_F = 0.527; T2_S = 0.0837; T2_F = 0.0166; M0_F = 0.23; k_FS = 6.7;

% GT-signals.
SPGR_Data = SPGR_steady_state_M0(FA_SPGR,TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'M0_F',M0_F,'k_FS',k_FS);
SSFP_Data = SSFP_steady_state_M0(FA_SSFP,TR_SSFP,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);
SSFP_Data_180 = SSFP_steady_state_180_M0(FA_SSFP,TR_SSFP,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);
Data = [SPGR_Data ; SSFP_Data; SSFP_Data_180];

Sigma = 1/300; % SNR = 300.

% Toby's Bounds.
Upper = [1.5 0.8 0.35 40 0.15 0.03];
Lower = [0.9 0.3 0.001 (1/0.6) 0.04 0.01];

%% Direct sampling of cost-function.

Steps = 200;

T1S_Vector = linspace(Lower(1),Upper(1),Steps);
T1F_Vector = linspace(Lower(2),Upper(2),Steps);
M0F_Vector = linspace(Lower(3),Upper(3),Steps);
kFS_Vector = linspace(Lower(4),Upper(4),Steps);
T2S_Vector = linspace(Lower(5),Upper(5),Steps);
T2F_Vector = linspace(Lower(6),Upper(6),Steps);

P_T1S = zeros(length(T1S_Vector),length(T1S_Vector));
P_T1F = zeros(length(T1S_Vector),length(T1S_Vector));
P_kFS = zeros(length(T1S_Vector),length(T1S_Vector));
P_T2S = zeros(length(T1S_Vector),length(T1S_Vector));
P_T2F = zeros(length(T1S_Vector),length(T1S_Vector));

disp (['Estimated Run Time: ' num2str(0.0081*Steps^2) ' seconds']);

for ii = 1:Steps
    tic
    for jj = 1:Steps
        
        P_T1S(ii,jj) = (logpdf_2P([T1S_Vector(ii) T1_F M0F_Vector(jj) k_FS T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        P_T1F(ii,jj) = (logpdf_2P([T1_S T1F_Vector(ii) M0F_Vector(jj) k_FS T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        P_kFS(ii,jj) = (logpdf_2P([T1_S T1_F M0F_Vector(ii) kFS_Vector(jj) T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        P_T2S(ii,jj) = (logpdf_2P([T1_S T1_F M0F_Vector(ii) k_FS T2S_Vector(jj) T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        P_T2F(ii,jj) = (logpdf_2P([T1_S T1_F M0F_Vector(ii) k_FS T2_S T2F_Vector(jj)], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        
    end
   toc 
end

Matrix_kFS = exp(P_kFS);
Matrix_T1F = exp(P_T1F); Matrix_T1S = exp(P_T1S);
Matrix_T2F = exp(P_T2F); Matrix_T2S = exp(P_T2S);
Threshold = 0.9999999;

[RowX_kFS, ColX_kFS] = find(Matrix_kFS > Threshold);
[RowX_T1F, ColX_T1F] = find(Matrix_T1F > Threshold);
[RowX_T1S, ColX_T1S] = find(Matrix_T1S > Threshold);
[RowX_T2F, ColX_T2F] = find(Matrix_T2F > Threshold);
[RowX_T2S, ColX_T2S] = find(Matrix_T2S > Threshold);

figure(5)
subplot(3,2,1)
imagesc([min(M0F_Vector) max(M0F_Vector)],[min(T1S_Vector) max(T1S_Vector)],exp(P_T1S)); colorbar; shading interp; xlabel('M_{0F}'); ylabel('T_{1S}'); hold on; plot(M0_F,T1_S,'w.', 'MarkerSize', 20);
hold on; scatter(M0F_Vector(ColX_T1S), T1S_Vector(RowX_T1S), 'rx')

subplot(3,2,2)
imagesc([min(M0F_Vector) max(M0F_Vector)],[min(T1F_Vector) max(T1F_Vector)],exp(P_T1F)); colorbar; shading interp; xlabel('M_{0F}'); ylabel('T_{1F}'); hold on; plot(M0_F,T1_F,'w.', 'MarkerSize', 20);
hold on; scatter(M0F_Vector(ColX_T1F), T1F_Vector(RowX_T1F), 'rx')

subplot(3,2,3)
imagesc([min(kFS_Vector) max(kFS_Vector)],[min(M0F_Vector) max(M0F_Vector)],exp(P_kFS)); colorbar; shading interp; xlabel('k_{FS}'); ylabel('M_{0F}'); hold on; plot(k_FS,M0_F,'w.', 'MarkerSize', 20);
hold on; scatter(kFS_Vector(ColX_kFS), M0F_Vector(RowX_kFS), 'rx')

subplot(3,2,4)
imagesc([min(T2S_Vector) max(T2S_Vector)],[min(M0F_Vector) max(M0F_Vector)],exp(P_T2S)); colorbar; shading interp; xlabel('T_{2S}'); ylabel('M_{0F}'); hold on; plot(T2_S,M0_F,'w.', 'MarkerSize', 20);
hold on; scatter(T2S_Vector(ColX_T2S), M0F_Vector(RowX_T2S), 'rx')

subplot(3,2,5)
imagesc([min(T2F_Vector) max(T2F_Vector)],[min(M0F_Vector) max(M0F_Vector)],exp(P_T2F)); colorbar; shading interp; xlabel('T_{2F}'); ylabel('M_{0F}'); hold on; plot(T2_F,M0_F,'w.', 'MarkerSize', 20);
hold on; scatter(T2F_Vector(ColX_T2F), M0F_Vector(RowX_T2F), 'rx')
