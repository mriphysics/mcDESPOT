%% Performs single-pool heat map simulations.

close all; clear all;

% Tissue and sequence parameters.
T1 = 1; T2 = 0.1; M0 = 1; Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;
TR_SPGR = 6.5e-3; TR_SSFP = 6.5e-3; FA_SPGR = deg2rad([2 4 6 8 10 12 14 16 18 20]); FA_SSFP180 = deg2rad([2 6 14 22 30 38 46 54 62 70]); FA_SSFP0 = deg2rad([2 6 14 22 30 38 46 54 62 70]);
%TR_SPGR = 7e-3; TR_SSFP = 7e-3; FA_SPGR = deg2rad([6 8 10 12 14 16]); FA_SSFP180 = deg2rad([15 25 35 45 55 65]); FA_SSFP0 = deg2rad([25 55]);

% Ground-truth signals for mcDESPOT.
SPGR_Data = SPGR_SP_SteadyState(FA_SPGR, TR_SPGR,'T1',T1,'M0',M0);
SSFP_Data_0 = SSFP_SP_SteadyState(FA_SSFP0, TR_SSFP, PC1,'T1',T1,'T2',T2,'M0',M0);
SSFP_Data_180 = SSFP_SP_SteadyState(FA_SSFP180, TR_SSFP, PC2,'T1',T1,'T2',T2,'M0',M0);

% Concatenate SSFP signals.
SSFP_Data = [SSFP_Data_0 ; SSFP_Data_180];

%% Direct sampling of cost-function in 2D. Modify for signal fit plots.
SNR = 30; Sigma = mean(SPGR_Data)/SNR;

Upper = [1.2 1.2 0.125]; Lower = [0.8 0.8 0.075];

% Pre-normalisation.
Data_O = [SPGR_Data ; SSFP_Data_0 ; SSFP_Data_180];

Steps = 50; nTrials = 50000;

T1_Vector = linspace(Lower(1),Upper(1),Steps); M0_Vector = linspace(Lower(2),Upper(2),Steps); 
T2_Vector = linspace(Lower(3),Upper(3),Steps); %Delta_Vector = linspace(Lower(4),Upper(4),Steps);

P_T1 = zeros(Steps,Steps,nTrials); P_T2 = zeros(Steps,Steps,nTrials); 
P_T1T2 = zeros(Steps,Steps,nTrials);

for ii = 1:Steps
    disp(['Outer Step Number: ', num2str(ii), '/', num2str(Steps), '.'])
    for jj = 1:Steps
        disp(['     Inner Step Number: ', num2str(jj), '/', num2str(Steps), '.'])
        tic
        parfor nn = 1:nTrials
            
            % Choose GT as first combination.
            %if nn == 1
                %T1_Rand = T1; T2_Rand = T2; M0_Rand = M0; Delta_Rand = Delta;
            %else
                
            % Draw random parameter values from a uniform distribution.
            T1_Rand = (Upper(1) - Lower(1)) .* rand(1,1) + Lower(1);
            M0_Rand = (Upper(2) - Lower(2)) .* rand(1,1) + Lower(2);
            T2_Rand = (Upper(3) - Lower(3)) .* rand(1,1) + Lower(3);
            %Delta_Rand = (Upper(4) - Lower(4)) .* rand(1,1) + Lower(4);
                
            %end
            
            P_T1(ii,jj,nn) = (logpdf_SinglePool([T1_Vector(ii) M0_Vector(jj) T2_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
            P_T2(ii,jj,nn) = (logpdf_SinglePool([T1_Rand M0_Vector(ii) T2_Vector(jj)], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
            P_T1T2(ii,jj,nn) = (logpdf_SinglePool([T1_Vector(ii) M0_Rand T2_Vector(jj)], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
            
        end
    end
end

%% Heat map analysis and plots.

figure(3); subplot(3,1,1);
MaxCF_T1 = max(exp(P_T1),[],3);
imagesc([min(M0_Vector) max(M0_Vector)],[min(T1_Vector) max(T1_Vector)],MaxCF_T1); colormap(viridis); colorbar; shading interp; xlabel('M_{0}','FontSize',12); ylabel('T_{1} (s)','FontSize',12); hold on;
plot(M0,T1,'w.', 'MarkerSize', 20);

subplot(3,1,2);
MaxCF_T2 = max(exp(P_T2),[],3);
imagesc([min(T2_Vector) max(T2_Vector)],[min(M0_Vector) max(M0_Vector)],MaxCF_T2); colormap(viridis); colorbar; shading interp; xlabel('T_{2} (s)','FontSize',12); ylabel('M_{0}','FontSize',12); hold on;
plot(T2,M0,'w.', 'MarkerSize', 20);

subplot(3,1,3)
MaxCF_T1T2 = max(exp(P_T1T2),[],3);
imagesc([min(T2_Vector) max(T2_Vector)],[min(T1_Vector) max(T1_Vector)],MaxCF_T1T2); colormap(viridis); colorbar; shading interp; xlabel('T_{2} (s)','FontSize',12); ylabel('T_{1} (s)','FontSize',12); hold on;
plot(T2,T1,'w.', 'MarkerSize', 20);
