clear all; close all

%% Choose GT signal set and concatenate.
 
TR_SPGR = 5e-3; TR_SSFP = 5e-3;
FA_SPGR = deg2rad([3 4 5 6 7 9 13 18]); FA_SSFP = deg2rad([10 13 17 20 23 30 43 60]);
T1_S = 0.965; T1_F = 0.465; T2_S = 0.09; T2_F = 0.012; M0_F = 0.2; k_FS = 8; % Deoni.

% GT-signals.
SPGR_Data = SPGR_SteadyState(FA_SPGR,TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'M0_F',M0_F,'k_FS',k_FS);
SSFP_Data = SSFP_SteadyState(FA_SSFP,TR_SSFP,0,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);
SSFP_Data_180 = SSFP_SteadyState(FA_SSFP,TR_SSFP,pi,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);
Data = [SPGR_Data ; SSFP_Data ; SSFP_Data_180];

%% Compute Hessian and perform eignevector decomposition.

Sig_SPGR = @(x)(SPGR_SteadyState(FA_SPGR,TR_SPGR,'T1_S',x(1),'T1_F',x(2),'M0_F',x(3),'k_FS',x(4)));
Sig_SSFP = @(x)(SSFP_SteadyState(FA_SSFP,TR_SSFP,0,'T1_S',x(1),'T2_S',x(5),'T1_F',x(2),'T2_F',x(6),'M0_F',x(3),'k_FS',x(4)));
Sig_SSFP180 = @(x)(SSFP_SteadyState(FA_SSFP,TR_SSFP,pi,'T1_S',x(1),'T2_S',x(5),'T1_F',x(2),'T2_F',x(6),'M0_F',x(3),'k_FS',x(4)));
Sig = @(x)[Sig_SPGR(x) ; Sig_SSFP(x) ; Sig_SSFP180(x)];

x0 = [T1_S T1_F M0_F k_FS T2_S T2_F];
CF = @(x)(norm(Sig(x) - Data).^2); 

% Treat x0 as the central parameter values, define scaling terms relative to these.
xs0 = ones(1,6);
CFs = @(xs)(norm(Sig(xs.*x0) - Data).^2);
Hessian = hessian(CFs, xs0);
[Eigenvectors,Eigenvalues] = eig(Hessian);

% Var = zeros(6,1);
% for nn = 1:6
%     Var(nn) = CF(x0 + Step * Eigenvectors(:,nn)');
% end

New_Coords = zeros(6,6); Scale = 1;
for pp = 1:6
    New_Coords(pp,:) = x0 + Scale * (Eigenvectors(:,pp).*x0')';
end

Coords_T1S = New_Coords(:,1); Coords_T1F = New_Coords(:,2); Coords_M0F = New_Coords(:,3); Coords_kFS = New_Coords(:,4); Coords_T2S = New_Coords(:,5); Coords_T2F = New_Coords(:,6);
Plot_Coords_M0FT1S = [[x0(3) ; Coords_M0F] , [x0(1) ; Coords_T1S]];
Plot_Coords_M0FT1F = [[x0(3) ; Coords_M0F] , [x0(2) ; Coords_T1F]];
Plot_Coords_M0FT2S = [[x0(5) ; Coords_T2S] , [x0(3) ; Coords_M0F]];
Plot_Coords_M0FT2F = [[x0(6) ; Coords_T2F] , [x0(3) ; Coords_M0F]];
Plot_Coords_M0FkFS = [[x0(4) ; Coords_kFS] , [x0(3) ; Coords_M0F]];
Plot_Coords_T1 = [[x0(2) ; Coords_T1F] , [x0(1) ; Coords_T1S]];
Plot_Coords_T2 = [[x0(6) ; Coords_T2F] , [x0(5) ; Coords_T2S]];
 
H_Inv = inv(Hessian);
Diag = diag(H_Inv); Diag_Scale = Diag ./ x0';

%% Direct sampling of cost-function.

Sigma = 1/300; % SNR = 30.
Upper = [2.0 1.0 0.500 50.00 0.200 0.050]; Lower = [0.5 0.1 0.001 0.001 0.001 0.001];
Steps = 200;
    
T1S_Vector = linspace(Lower(1),Upper(1),Steps); T1F_Vector = linspace(Lower(2),Upper(2),Steps);
M0F_Vector = linspace(Lower(3),Upper(3),Steps); kFS_Vector = linspace(Lower(4),Upper(4),Steps);
T2S_Vector = linspace(Lower(5),Upper(5),Steps); T2F_Vector = linspace(Lower(6),Upper(6),Steps);

P_T1S = zeros(length(M0_F),length(T1S_Vector),length(T1S_Vector));
P_T1F = zeros(length(M0_F),length(T1S_Vector),length(T1S_Vector));
P_kFS = zeros(length(M0_F),length(T1S_Vector),length(T1S_Vector));
P_T2S = zeros(length(M0_F),length(T1S_Vector),length(T1S_Vector));
P_T2F = zeros(length(M0_F),length(T1S_Vector),length(T1S_Vector));
P_T1 = zeros(length(M0_F),length(T1S_Vector),length(T1S_Vector));
P_T2 = zeros(length(M0_F),length(T1S_Vector),length(T1S_Vector));

for ii = 1:Steps
    tic
    for jj = 1:Steps
        
        % Uses results of eigenvector decomposition to define cut.
        %P_T1S(ii,jj) = (logpdf([T1S_Vector(ii) Coords_T1F(6) M0F_Vector(jj) Coords_kFS(6) Coords_T2S(6) Coords_T2F(6)], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        %P_T1F(ii,jj) = (logpdf([Coords_T1S(6) T1F_Vector(ii) M0F_Vector(jj) Coords_kFS(6) Coords_T2S(6) Coords_T2F(6)], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        %P_kFS(ii,jj) = (logpdf([Coords_T1S(6) Coords_T1F(6) M0F_Vector(ii) kFS_Vector(jj) Coords_T2S(6) Coords_T2F(6)], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        %P_T2S(ii,jj) = (logpdf([Coords_T1S(6) Coords_T1F(6) M0F_Vector(ii) Coords_kFS(6) T2S_Vector(jj) Coords_T2F(6)], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        %P_T2F(ii,jj) = (logpdf([Coords_T1S(6) Coords_T1F(6) M0F_Vector(ii) Coords_kFS(6) Coords_T2S(6) T2F_Vector(jj)], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        %P_T1(ii,jj)  = (logpdf([T1S_Vector(ii) T1F_Vector(jj) Coords_M0F(6) Coords_kFS(6) Coords_T2S(6) Coords_T2F(6)], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        %P_T2(ii,jj)  = (logpdf([Coords_T1S(6) Coords_T1F(6) Coords_M0F(6) Coords_kFS(6) T2S_Vector(ii) T2F_Vector(jj)], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        
        % Uses GT to define cut.
        P_T1S(ii,jj) = (logpdf([T1S_Vector(ii) T1_F M0F_Vector(jj) k_FS T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        P_T1F(ii,jj) = (logpdf([T1_S T1F_Vector(ii) M0F_Vector(jj) k_FS T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        P_kFS(ii,jj) = (logpdf([T1_S T1_F M0F_Vector(ii) kFS_Vector(jj) T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        P_T2S(ii,jj) = (logpdf([T1_S T1_F M0F_Vector(ii) k_FS T2S_Vector(jj) T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        P_T2F(ii,jj) = (logpdf([T1_S T1_F M0F_Vector(ii) k_FS T2_S T2F_Vector(jj)], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        P_T1(ii,jj)  = (logpdf([T1S_Vector(ii) T1F_Vector(jj) M0_F k_FS T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        P_T2(ii,jj)  = (logpdf([T1_S T1_F M0_F k_FS T2S_Vector(ii) T2F_Vector(jj)], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        
    end
    toc
end

%% Perform fitting using stochastic region contraction.

Trials = 5000; Iterations = 7; N = 50; Runs = 1; Params = 6;
[Solution, Bounds, Subsets] = SRC_Sim(Params, Trials, Iterations, N, Runs, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data);

% Bounds Order: T1S_LB, T1S_UB, T1F_LB, T1F_UB, T2S_LB, T2S_UB, T2F_LB, T2F_UB, M0F_LB, M0F_UB, kFS_LB, kFS_UB];
T1S_LB = Bounds(:,1); T1S_UB = Bounds(:,2); T1F_LB = Bounds(:,3); T1F_UB = Bounds(:,4);
T2S_LB = Bounds(:,5); T2S_UB = Bounds(:,6); T2F_LB = Bounds(:,7); T2F_UB = Bounds(:,8);
M0F_LB = Bounds(:,9); M0F_UB = Bounds(:,10); kFS_LB = Bounds(:,11); kFS_UB = Bounds(:,12);

%% Plot results: fake_parula, viridis, plasma, magma, inferno.

LineSettings = ['w', 'w', 'w', 'w', 'w', 'w', 'w', 'w'];
ArrowColor = ['w','b','c','g','y','m','r'];

figure(5)
subplot(3,2,1)
imagesc([min(M0F_Vector) max(M0F_Vector)],[min(T1S_Vector) max(T1S_Vector)],exp(squeeze(P_T1S))); colormap(magma); colorbar; shading interp; xlabel('M_{0F}'); ylabel('T_{1S}'); 
hold on;
% for aa = 2:7
%     arrow(Plot_Coords_M0FT1S(1,:),Plot_Coords_M0FT1S(aa,:),'Length',1,'Color', ArrowColor(aa));
% end
% ll = legend('V_{1}','V_{2}','V_{3}','V_{4}','V_{5}','V_{6}'); ll.FontSize = 12;
plot(M0_F,T1_S,'w.', 'MarkerSize', 30);
 
%arrow(Plot_Coords_M0FT1S(1,:),Plot_Coords_M0FT1S(3,:),'Length',1,'Color', 'b')
%arrow(Plot_Coords_M0FT1S(1,:),Plot_Coords_M0FT1S(6,:),'Length',1,'Color', 'r')

% for rr = 1:length(T1S_LB)
%     rectangle('Position',[M0F_LB(rr) T1S_LB(rr) (M0F_UB(rr)-M0F_LB(rr)) (T1S_UB(rr)-T1S_LB(rr))],'EdgeColor',LineSettings(rr)) 
% end
% for ss = 1:size(Subsets,1)
%        plot(Subsets(ss,5),Subsets(ss,1),'r+')
% end

subplot(3,2,2)
imagesc([min(M0F_Vector) max(M0F_Vector)],[min(T1F_Vector) max(T1F_Vector)],exp(squeeze(P_T1F))); colormap(magma); colorbar; shading interp; xlabel('M_{0F}'); ylabel('T_{1F}'); 
hold on; 
% for aa = 2:7
%     arrow(Plot_Coords_M0FT1F(1,:),Plot_Coords_M0FT1F(aa,:),'Length',1,'Color', ArrowColor(aa))
% end
% ll = legend('V_{1}','V_{2}','V_{3}','V_{4}','V_{5}','V_{6}'); ll.FontSize = 12;
plot(M0_F,T1_F,'w.', 'MarkerSize',30);

%arrow(Plot_Coords_M0FT1F(1,:),Plot_Coords_M0FT1F(2,:),'Length',1,'Color', 'b')
%arrow(Plot_Coords_M0FT1F(1,:),Plot_Coords_M0FT1F(7,:),'Length',1,'Color', 'r')

% for rr = 1:length(T1S_LB)
%     rectangle('Position',[M0F_LB(rr) T1F_LB(rr) (M0F_UB(rr)-M0F_LB(rr)) (T1F_UB(rr)-T1F_LB(rr))],'EdgeColor',LineSettings(rr)) 
% end
% for ss = 1:size(Subsets,1)
%        plot(Subsets(ss,5),Subsets(ss,2),'r+')
% end

subplot(3,2,3)
imagesc([min(kFS_Vector) max(kFS_Vector)],[min(M0F_Vector) max(M0F_Vector)],exp(squeeze(P_kFS))); colormap(magma); colorbar; shading interp; xlabel('k_{FS}'); ylabel('M_{0F}'); 
hold on; 
% for aa = 2:7
%     arrow(Plot_Coords_M0FkFS(1,:),Plot_Coords_M0FkFS(aa,:),'Length',1,'Color', ArrowColor(aa))
% end
% ll = legend('V_{1}','V_{2}','V_{3}','V_{4}','V_{5}','V_{6}'); ll.FontSize = 12;
plot(k_FS,M0_F,'w.', 'MarkerSize', 30);

%arrow(Plot_Coords_M0FkFS(1,:),Plot_Coords_M0FkFS(2,:),'Length',1,'Color', 'b')
%arrow(Plot_Coords_M0FkFS(1,:),Plot_Coords_M0FkFS(7,:),'Length',1,'Color', 'r')

% for rr = 1:length(T1S_LB)
%     rectangle('Position',[kFS_LB(rr) M0F_LB(rr) (kFS_UB(rr)-kFS_LB(rr)) (M0F_UB(rr)-M0F_LB(rr))],'EdgeColor',LineSettings(rr)) 
% end
% for ss = 1:size(Subsets,1)
%        plot(Subsets(ss,6),Subsets(ss,5),'r+')
% end

subplot(3,2,4)
imagesc([min(T2S_Vector) max(T2S_Vector)],[min(M0F_Vector) max(M0F_Vector)],exp(squeeze(P_T2S))); colormap(magma); colorbar; shading interp; xlabel('T_{2S}'); ylabel('M_{0F}'); 
hold on; 
% for aa = 2:7
%     arrow(Plot_Coords_M0FT2S(1,:),Plot_Coords_M0FT2S(aa,:),'Length',1,'Color', ArrowColor(aa))
% end
% ll = legend('V_{1}','V_{2}','V_{3}','V_{4}','V_{5}','V_{6}'); ll.FontSize = 12;
plot(T2_S,M0_F,'w.', 'MarkerSize', 30);

% arrow(Plot_Coords_M0FT2S(1,:),Plot_Coords_M0FT2S(2,:),'Length',1,'Color', 'b')
% arrow(Plot_Coords_M0FT2S(1,:),Plot_Coords_M0FT2S(7,:),'Length',1,'Color', 'r')

% for rr = 1:length(T1S_LB)
%     rectangle('Position',[T2S_LB(rr) M0F_LB(rr) (T2S_UB(rr)-T2S_LB(rr)) (M0F_UB(rr)-M0F_LB(rr))],'EdgeColor',LineSettings(rr)) 
% end
% for ss = 1:size(Subsets,1)
%        plot(Subsets(ss,3),Subsets(ss,5),'r+')
% end

subplot(3,2,5)
imagesc([min(T2F_Vector) max(T2F_Vector)],[min(M0F_Vector) max(M0F_Vector)],exp(squeeze(P_T2F))); colormap(magma); colorbar; shading interp; xlabel('T_{2F}'); ylabel('M_{0F}'); 
hold on; 
% for aa = 2:7
%     arrow(Plot_Coords_M0FT2F(1,:),Plot_Coords_M0FT2F(aa,:),'Length',1,'Color', ArrowColor(aa))
% end
% ll = legend('V_{1}','V_{2}','V_{3}','V_{4}','V_{5}','V_{6}'); ll.FontSize = 12;
plot(T2_F,M0_F,'w.', 'MarkerSize', 30);

%arrow(Plot_Coords_M0FT2F(1,:),Plot_Coords_M0FT2F(2,:),'Length',1,'Color', 'b')
%arrow(Plot_Coords_M0FT2F(1,:),Plot_Coords_M0FT2F(7,:),'Length',1,'Color', 'r')

% for rr = 1:length(T1S_LB)
%     rectangle('Position',[T2F_LB(rr) M0F_LB(rr) (T2F_UB(rr)-T2F_LB(rr)) (M0F_UB(rr)-M0F_LB(rr))],'EdgeColor',LineSettings(rr)) 
% end
% for ss = 1:size(Subsets,1)
%        plot(Subsets(ss,4),Subsets(ss,5),'r+')
% end

% imagesc([min(T1F_Vector) max(T1F_Vector)],[min(T1S_Vector) max(T1S_Vector)],exp(squeeze(P_T1))); colormap(magma); colorbar; shading interp; xlabel('T_{1F}'); ylabel('T_{1S}'); 
% hold on;
% for aa = 2:7
%     arrow(Plot_Coords_T1(1,:),Plot_Coords_T1(aa,:),'Length',1,'Color', ArrowColor(aa));
% end
% ll = legend('V_{1}','V_{2}','V_{3}','V_{4}','V_{5}','V_{6}'); ll.FontSize = 12;
% plot(T1_F,T1_S,'w.', 'MarkerSize', 30);
% 
% imagesc([min(T2F_Vector) max(T2F_Vector)],[min(T2S_Vector) max(T2S_Vector)],exp(squeeze(P_T2))); colormap(magma); colorbar; shading interp; xlabel('T_{2F}'); ylabel('T_{2S}'); 
% hold on;
% for aa = 2:7
%     arrow(Plot_Coords_T2(1,:),Plot_Coords_T2(aa,:),'Length',1,'Color', ArrowColor(aa));
% end
% ll = legend('V_{1}','V_{2}','V_{3}','V_{4}','V_{5}','V_{6}'); ll.FontSize = 12;
% plot(T2_F,T2_S,'w.', 'MarkerSize', 30);