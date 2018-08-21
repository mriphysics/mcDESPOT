TR_SPGR = 6.5e-3; TR_SSFP = 6.5e-3; FA_SPGR = deg2rad([2 4 6 8 10 12 14 16 18 20]); FA_SSFP180 = deg2rad([2 6 14 22 30 38 46 54 62 70]); FA_SSFP0 = deg2rad([2 6 14 22 30 38 46 54 62 70]);
T1S_1 = 1; T1F_1 = 0.35; T2S_1 = 0.080; T2F_1 = 0.015; kFS_1 = 10; M0F_1 = 0.25; Delta_1 = 0; PC1_1 = 0 + Delta_1; PC2_1 = pi + Delta_1;
T1S_2 = 1.15; T1F_2 = 0.4; T2S_2 = 0.11; T2F_2 = 0.02; kFS_2 = 7.5; M0F_2 = 0.175; Delta_2 = 0; PC1_2 = 0 + Delta_2; PC2_2 = pi + Delta_2;
T1S_3 = 1.3; T1F_3 = 0.45; T2S_3 = 0.14; T2F_3 = 0.025; kFS_3 = 5; M0F_3 = 0.1; Delta_3 = 0; PC1_3 = 0 + Delta_3; PC2_3 = pi + Delta_3;

%T1S_1 = 1; T1F_1 = 0.35; T2S_1 = 0.080; T2F_1 = 0.015; kFS_1 = 0; M0F_1 = 0.25; Delta_1 = 0; PC1_1 = 0 + Delta_1; PC2_1 = pi + Delta_1;
%T1S_2 = 1.15; T1F_2 = 0.4; T2S_2 = 0.11; T2F_2 = 0.02; kFS_2 = 0; M0F_2 = 0.175; Delta_2 = 0; PC1_2 = 0 + Delta_2; PC2_2 = pi + Delta_2;
%T1S_3 = 1.3; T1F_3 = 0.45; T2S_3 = 0.14; T2F_3 = 0.025; kFS_3 = 0; M0F_3 = 0.1; Delta_3 = 0; PC1_3 = 0 + Delta_3; PC2_3 = pi + Delta_3;

%T1S_1 = 1; T1F_1 = 0.35; T2S_1 = 0.080; T2F_1 = 0.015; kFS_1 = 5; M0F_1 = 0.25; Delta_1 = 0; PC1_1 = 0 + Delta_1; PC2_1 = pi + Delta_1;
%T1S_2 = 1; T1F_2 = 0.35; T2S_2 = 0.080; T2F_3 = 0.015; kFS_2 = 7.5; M0F_2 = 0.25; Delta_2 = 0; PC1_2 = 0 + Delta_2; PC2_2 = pi + Delta_2;
%T1S_3 = 1; T1F_3 = 0.35; T2S_3 = 0.080; T2F_2 = 0.015; kFS_3 = 10; M0F_3 = 0.25; Delta_3 = 0; PC1_3 = 0 + Delta_3; PC2_3 = pi + Delta_3;

nTrials = 200e6; PlottingNo = 1000;

SPGR_Data_1 = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1S_1,'T1_F',T1F_1,'M0_F',M0F_1,'k_FS',kFS_1);
SSFP_Data_0_1 = SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1_1,'T1_S',T1S_1,'T2_S',T2S_1,'T1_F',T1F_1,'T2_F',T2F_1,'M0_F',M0F_1,'k_FS',kFS_1);
SSFP_Data_180_1 = SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2_1,'T1_S',T1S_1,'T2_S',T2S_1,'T1_F',T1F_1,'T2_F',T2F_1,'M0_F',M0F_1,'k_FS',kFS_1);
SSFP_Data_1 = [SSFP_Data_0_1 ; SSFP_Data_180_1];

SPGR_Data_2 = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1S_2,'T1_F',T1F_2,'M0_F',M0F_2,'k_FS',kFS_2);
SSFP_Data_0_2 = SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1_2,'T1_S',T1S_2,'T2_S',T2S_2,'T1_F',T1F_2,'T2_F',T2F_2,'M0_F',M0F_2,'k_FS',kFS_2);
SSFP_Data_180_2 = SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2_2,'T1_S',T1S_2,'T2_S',T2S_2,'T1_F',T1F_2,'T2_F',T2F_2,'M0_F',M0F_2,'k_FS',kFS_2);
SSFP_Data_2 = [SSFP_Data_0_2 ; SSFP_Data_180_2];

SPGR_Data_3 = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1S_3,'T1_F',T1F_3,'M0_F',M0F_3,'k_FS',kFS_3);
SSFP_Data_0_3 = SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1_3,'T1_S',T1S_3,'T2_S',T2S_3,'T1_F',T1F_3,'T2_F',T2F_3,'M0_F',M0F_3,'k_FS',kFS_3);
SSFP_Data_180_3 = SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2_3,'T1_S',T1S_3,'T2_S',T2S_3,'T1_F',T1F_3,'T2_F',T2F_3,'M0_F',M0F_3,'k_FS',kFS_3);
SSFP_Data_3 = [SSFP_Data_0_3 ; SSFP_Data_180_3];

figure(1);
plot(rad2deg(FA_SPGR), SPGR_Data_1/mean(SPGR_Data_1), 'bo','Linewidth',2,'MarkerSize',10); hold on
plot(rad2deg(FA_SSFP0), SSFP_Data_0_1/mean(SSFP_Data_1), 'ko','LineWidth',2,'MarkerSize',10)
plot(rad2deg(FA_SSFP180), SSFP_Data_180_1/mean(SSFP_Data_1), 'ro','LineWidth',2,'MarkerSize',10)
plot(rad2deg(FA_SPGR), SPGR_Data_2/mean(SPGR_Data_2), 'bo','Linewidth',2,'MarkerSize',10);
plot(rad2deg(FA_SSFP0), SSFP_Data_0_2/mean(SSFP_Data_2), 'ko','LineWidth',2,'MarkerSize',10)
plot(rad2deg(FA_SSFP180), SSFP_Data_180_2/mean(SSFP_Data_2), 'ro','LineWidth',2,'MarkerSize',10)
plot(rad2deg(FA_SPGR), SPGR_Data_3/mean(SPGR_Data_3), 'bo','Linewidth',2,'MarkerSize',10);
plot(rad2deg(FA_SSFP0), SSFP_Data_0_3/mean(SSFP_Data_3), 'ko','LineWidth',2,'MarkerSize',10)
plot(rad2deg(FA_SSFP180), SSFP_Data_180_3/mean(SSFP_Data_3), 'ro','LineWidth',2,'MarkerSize',10)
xlabel('FA [^{o}]', 'FontSize', 14); ylabel('Normalised Signal (a.u.)', 'FontSize',14);
ll = legend('SPGR','bSSFP_{0}','bSSFP_{180}'); ll.FontSize = 16; ll.AutoUpdate = 'off'; legend('boxoff'); grid on; grid minor;
 
SPGR_FM_1 = zeros(length(FA_SPGR),PlottingNo);
SSFP0_FM_1 = zeros(length(FA_SSFP0),PlottingNo);
SSFP180_FM_1 = zeros(length(FA_SSFP180),PlottingNo);
SSFP_FM_1 = zeros((length(FA_SSFP0)+length(FA_SSFP180)),PlottingNo);
SPGR_FM_2 = zeros(length(FA_SPGR),PlottingNo);
SSFP0_FM_2 = zeros(length(FA_SSFP0),PlottingNo);
SSFP180_FM_2 = zeros(length(FA_SSFP180),PlottingNo);
SSFP_FM_2 = zeros((length(FA_SSFP0)+length(FA_SSFP180)),PlottingNo);
SPGR_FM_3 = zeros(length(FA_SPGR),PlottingNo);
SSFP0_FM_3 = zeros(length(FA_SSFP0),PlottingNo);
SSFP180_FM_3 = zeros(length(FA_SSFP180),PlottingNo);
SSFP_FM_3 = zeros((length(FA_SSFP0)+length(FA_SSFP180)),PlottingNo);

for ii = 1:PlottingNo

    SPGR_FM_1(:,ii) = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1S_Picked_1(ii),'T1_F',T1F_Picked_1(ii),'M0_F',M0F_Picked_1(ii),'k_FS',0);
    SSFP0_FM_1(:,ii) = SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1_1,'T1_S',T1S_Picked_1(ii),'T2_S',T2S_Picked_1(ii),'T1_F',T1F_Picked_1(ii),'T2_F',T2F_Picked_1(ii),'M0_F',M0F_Picked_1(ii),'k_FS',0);
    SSFP180_FM_1(:,ii) = SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2_1,'T1_S',T1S_Picked_1(ii),'T2_S',T2S_Picked_1(ii),'T1_F',T1F_Picked_1(ii),'T2_F',T2F_Picked_1(ii),'M0_F',M0F_Picked_1(ii),'k_FS',0);
    SSFP_FM_1(:,ii) = [SSFP0_FM_1(:,ii) ; SSFP180_FM_1(:,ii)];
    
    SPGR_FM_2(:,ii) = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1S_Picked_2(ii),'T1_F',T1F_Picked_2(ii),'M0_F',M0F_Picked_2(ii),'k_FS',0);
    SSFP0_FM_2(:,ii) = SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1_2,'T1_S',T1S_Picked_2(ii),'T2_S',T2S_Picked_2(ii),'T1_F',T1F_Picked_2(ii),'T2_F',T2F_Picked_2(ii),'M0_F',M0F_Picked_2(ii),'k_FS',0);
    SSFP180_FM_2(:,ii) = SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2_2,'T1_S',T1S_Picked_2(ii),'T2_S',T2S_Picked_2(ii),'T1_F',T1F_Picked_2(ii),'T2_F',T2F_Picked_2(ii),'M0_F',M0F_Picked_2(ii),'k_FS',0);
    SSFP_FM_2(:,ii) = [SSFP0_FM_2(:,ii) ; SSFP180_FM_2(:,ii)];
    
    SPGR_FM_3(:,ii) = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1S_Picked_3(ii),'T1_F',T1F_Picked_3(ii),'M0_F',M0F_Picked_3(ii),'k_FS',0);
    SSFP0_FM_3(:,ii) = SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1_3,'T1_S',T1S_Picked_3(ii),'T2_S',T2S_Picked_3(ii),'T1_F',T1F_Picked_3(ii),'T2_F',T2F_Picked_3(ii),'M0_F',M0F_Picked_3(ii),'k_FS',0);
    SSFP180_FM_3(:,ii) = SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2_3,'T1_S',T1S_Picked_3(ii),'T2_S',T2S_Picked_3(ii),'T1_F',T1F_Picked_3(ii),'T2_F',T2F_Picked_3(ii),'M0_F',M0F_Picked_3(ii),'k_FS',0);
    SSFP_FM_3(:,ii) = [SSFP0_FM_3(:,ii) ; SSFP180_FM_3(:,ii)];
    
end
 

for jj = 1:PlottingNo

    figure(1); plot(rad2deg(FA_SPGR), SPGR_FM_1(:,jj)/mean(SPGR_FM_1(:,jj)), '--', 'Color', 'c', 'LineWidth', 0.1);
    plot(rad2deg(FA_SSFP0), SSFP0_FM_1(:,jj)/mean(SSFP_FM_1(:,jj)), '--', 'Color', 'c', 'LineWidth',0.1);
    plot(rad2deg(FA_SSFP180), SSFP180_FM_1(:,jj)/mean(SSFP_FM_1(:,jj)), '--', 'Color', 'c', 'LineWidth',0.1); hold on;
    
    figure(1); plot(rad2deg(FA_SPGR), SPGR_FM_2(:,jj)/mean(SPGR_FM_2(:,jj)), '--', 'Color', 'k', 'LineWidth', 0.1);
    plot(rad2deg(FA_SSFP0), SSFP0_FM_2(:,jj)/mean(SSFP_FM_2(:,jj)), '--', 'Color', 'k', 'LineWidth',0.1);
    plot(rad2deg(FA_SSFP180), SSFP180_FM_2(:,jj)/mean(SSFP_FM_2(:,jj)), '--', 'Color', 'k', 'LineWidth',0.1); hold on;
    
    figure(1); plot(rad2deg(FA_SPGR), SPGR_FM_3(:,jj)/mean(SPGR_FM_3(:,jj)), '--', 'Color', 'm', 'LineWidth', 0.1);
    plot(rad2deg(FA_SSFP0), SSFP0_FM_3(:,jj)/mean(SSFP_FM_3(:,jj)), '--', 'Color', 'm', 'LineWidth',0.1);
    plot(rad2deg(FA_SSFP180), SSFP180_FM_3(:,jj)/mean(SSFP_FM_3(:,jj)), '--', 'Color', 'm', 'LineWidth',0.1); hold on;

end

figure(2);
subplot(3,2,1); histogram(T1F_Picked_1,'BinWidth',0.01,'FaceColor','c','FaceAlpha',0.7); xlabel('T_{1F} (s)','FontSize',12); ylabel('Count','FontSize',12); hold on; 
histogram(T1F_Picked_2,'BinWidth',0.01,'FaceColor','k','FaceAlpha',0.7); 
histogram(T1F_Picked_3,'BinWidth',0.01,'FaceColor','m','FaceAlpha',0.7); 
ll = legend('Tissue 1','Tissue 2','Tissue 3'); ll.FontSize = 14; ll.AutoUpdate = 'off'; ll.Location = 'northwest';
line([T1F_1 T1F_1],[0 300],'Color','r','LineStyle','--','LineWidth',2); line([T1F_2 T1F_2],[0 300],'Color','r','LineStyle','--','LineWidth',2);  line([T1F_3 T1F_3],[0 300],'Color','r','LineStyle','--','LineWidth',2)

subplot(3,2,2); histogram(T1S_Picked_1,'BinWidth',0.005,'FaceColor','c','FaceAlpha',0.7); xlabel('T_{1S} (s)','FontSize',12); ylabel('Count','FontSize',12); hold on; line([T1S_1 T1S_1],[0 300],'Color','r','LineStyle','--','LineWidth',2)
histogram(T1S_Picked_2,'BinWidth',0.005,'FaceColor','k','FaceAlpha',0.7); line([T1S_2 T1S_2],[0 300],'Color','r','LineStyle','--','LineWidth',2)
histogram(T1S_Picked_3,'BinWidth',0.005,'FaceColor','m','FaceAlpha',0.7); line([T1S_3 T1S_3],[0 300],'Color','r','LineStyle','--','LineWidth',2)

subplot(3,2,3); histogram(T2F_Picked_1,'BinWidth',0.0003,'FaceColor','c','FaceAlpha',0.7); xlabel('T_{2F} (s)','FontSize',12); ylabel('Count','FontSize',12); hold on; line([T2F_1 T2F_1],[0 300],'Color','r','LineStyle','--','LineWidth',2)
histogram(T2F_Picked_2,'BinWidth',0.0003,'FaceColor','k','FaceAlpha',0.7); line([T2F_2 T2F_2],[0 300],'Color','r','LineStyle','--','LineWidth',2)
histogram(T2F_Picked_3,'BinWidth',0.0003,'FaceColor','m','FaceAlpha',0.7); line([T2F_3 T2F_3],[0 300],'Color','r','LineStyle','--','LineWidth',2)

subplot(3,2,4); histogram(T2S_Picked_1,'BinWidth',0.0005,'FaceColor','c','FaceAlpha',0.7); xlabel('T_{2S} (s)','FontSize',12); ylabel('Count','FontSize',12); hold on; line([T2S_1 T2S_1],[0 200],'Color','r','LineStyle','--','LineWidth',2)
histogram(T2S_Picked_2,'BinWidth',0.0005,'FaceColor','k','FaceAlpha',0.7); line([T2S_2 T2S_2],[0 200],'Color','r','LineStyle','--','LineWidth',2)
histogram(T2S_Picked_3,'BinWidth',0.0005,'FaceColor','m','FaceAlpha',0.7); line([T2S_3 T2S_3],[0 200],'Color','r','LineStyle','--','LineWidth',2)

subplot(3,2,5); histogram(M0F_Picked_1,'BinWidth',0.003,'FaceColor','c','FaceAlpha',0.7); xlabel('M_{0F}','FontSize',12); ylabel('Count','FontSize',12); hold on; line([M0F_1 M0F_1],[0 200],'Color','r','LineStyle','--','LineWidth',2)
histogram(M0F_Picked_2,'BinWidth',0.003,'FaceColor','k','FaceAlpha',0.7); line([M0F_2 M0F_2],[0 200],'Color','r','LineStyle','--','LineWidth',2)
histogram(M0F_Picked_3,'BinWidth',0.003,'FaceColor','m','FaceAlpha',0.7); line([M0F_3 M0F_3],[0 200],'Color','r','LineStyle','--','LineWidth',2)

% subplot(3,2,6); histogram(kFS_Picked_1,'BinWidth',2,'FaceColor','c','FaceAlpha',0.7); xlabel('k_{FS} (s^{-1})','FontSize',12); ylabel('Count','FontSize',12); hold on; line([kFS_1 kFS_1],[0 200],'Color','r','LineStyle','--','LineWidth',2)
% histogram(kFS_Picked_2,'BinWidth',2,'FaceColor','k','FaceAlpha',0.7); line([kFS_2 kFS_2],[0 200],'Color','r','LineStyle','--','LineWidth',2)
% histogram(kFS_Picked_3,'BinWidth',2,'FaceColor','m','FaceAlpha',0.7); line([kFS_3 kFS_3],[0 200],'Color','r','LineStyle','--','LineWidth',2)
