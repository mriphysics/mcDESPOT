%%% Generates signal plots and histograms for Figure 1. %%

TR_SPGR = 6.5e-3; TR_SSFP = 6.5e-3; FA_SPGR = deg2rad([2 4 6 8 10 12 14 16 18 20]); FA_SSFP180 = deg2rad([2 6 14 22 30 38 46 54 62 70]); FA_SSFP0 = deg2rad([2 6 14 22 30 38 46 54 62 70]);

exchange = 'on';

switch exchange
    case 'on'
        % HB.
        T1F_0 = 0.45; T1S_0 = 1.4; T2F_0 = 0.015; T2S_0 = 0.09; M0F_0 = 0.15; kFS_0 = 8; Delta_0 = 0; PC1_0 = 0 + Delta_0; PC2_0 = pi + Delta_0;
        % WML.
        T1S_1 = 1; T1F_1 = 0.35; T2S_1 = 0.080; T2F_1 = 0.015; kFS_1 = 10; M0F_1 = 0.25; Delta_1 = 0; PC1_1 = 0 + Delta_1; PC2_1 = pi + Delta_1;
        % INT.
        T1S_2 = 1.15; T1F_2 = 0.4; T2S_2 = 0.11; T2F_2 = 0.02; kFS_2 = 7.5; M0F_2 = 0.175; Delta_2 = 0; PC1_2 = 0 + Delta_2; PC2_2 = pi + Delta_2;
        % GML.
        T1S_3 = 1.3; T1F_3 = 0.45; T2S_3 = 0.14; T2F_3 = 0.025; kFS_3 = 5; M0F_3 = 0.1; Delta_3 = 0; PC1_3 = 0 + Delta_3; PC2_3 = pi + Delta_3;
    case 'off'
        % HB.
        T1F_0 = 0.45; T1S_0 = 1.4; T2F_0 = 0.015; T2S_0 = 0.09; M0F_0 = 0.15; kFS_0 = 0; Delta_0 = 0; PC1_0 = 0 + Delta_0; PC2_0 = pi + Delta_0;
        % WML.
        T1S_1 = 1; T1F_1 = 0.35; T2S_1 = 0.080; T2F_1 = 0.015; kFS_1 = 0; M0F_1 = 0.25; Delta_1 = 0; PC1_1 = 0 + Delta_1; PC2_1 = pi + Delta_1;
        % INT.
        T1S_2 = 1.15; T1F_2 = 0.4; T2S_2 = 0.11; T2F_2 = 0.02; kFS_2 = 0; M0F_2 = 0.175; Delta_2 = 0; PC1_2 = 0 + Delta_2; PC2_2 = pi + Delta_2;
        % GML.
        T1S_3 = 1.3; T1F_3 = 0.45; T2S_3 = 0.14; T2F_3 = 0.025; kFS_3 = 0; M0F_3 = 0.1; Delta_3 = 0; PC1_3 = 0 + Delta_3; PC2_3 = pi + Delta_3;
end
        
nTrials = 200e6; PlottingNo = 1000;

%% Histograms.

figure(1);
subplot(3,2,1); 
histogram(T1F_Picked_T0,'BinWidth',0.02,'DisplayStyle','stairs','LineWidth',2); hold on; xlabel('T_{1F} (s)','FontSize',18); ylabel('Count','FontSize',18);
histogram(T1F_Picked_T1,'BinWidth',0.02,'DisplayStyle','stairs','LineWidth',2); 
histogram(T1F_Picked_T2,'BinWidth',0.02,'DisplayStyle','stairs','LineWidth',2); 
histogram(T1F_Picked_T3,'BinWidth',0.02,'DisplayStyle','stairs','LineWidth',2); 
cm = get(gca,'ColorOrder');
line([T1F_1 T1F_1],[0 150],'Color',cm(2,:),'LineStyle',':','LineWidth',4); line([T1F_2 T1F_2],[0 150],'Color',cm(3,:),'LineStyle',':','LineWidth',4);  line([T1F_3 T1F_3],[0 150],'Color',cm(4,:),'LineStyle',':','LineWidth',4); line([T1F_0 T1F_0],[0 150],'Color',cm(1,:),'LineStyle',':','LineWidth',4)
get(gca, 'XTick'); set(gca, 'FontSize', 18); get(gca, 'YTick'); set(gca, 'FontSize', 18); ylim([0 150])
grid on; grid minor

subplot(3,2,2); 
histogram(T1S_Picked_T0,'BinWidth',0.05,'DisplayStyle','stairs','LineWidth',2); hold on 
histogram(T1S_Picked_T1,'BinWidth',0.05,'DisplayStyle','stairs','LineWidth',2); hold on; xlabel('T_{1S} (s)','FontSize',18); ylabel('Count','FontSize',18); 
histogram(T1S_Picked_T2,'BinWidth',0.05,'DisplayStyle','stairs','LineWidth',2); 
histogram(T1S_Picked_T3,'BinWidth',0.05,'DisplayStyle','stairs','LineWidth',2); 
ll = legend({'HB','WML','INT','GML'},'Position',[0.846265938069217,0.76394844802434,0.056284152158622,0.191433561014963]); ll.FontSize = 20; ll.AutoUpdate = 'off'; legend('boxoff');
line([T1S_1 T1S_1],[0 320],'Color',cm(2,:),'LineStyle',':','LineWidth',4); line([T1S_2 T1S_2],[0 320],'Color',cm(3,:),'LineStyle',':','LineWidth',4);  line([T1S_3 T1S_3],[0 320],'Color',cm(4,:),'LineStyle',':','LineWidth',4); line([T1S_0 T1S_0],[0 320],'Color',cm(1,:),'LineStyle',':','LineWidth',4)
get(gca, 'XTick'); set(gca, 'FontSize', 18); get(gca, 'YTick'); set(gca, 'FontSize', 18); ylim([0 320]);
grid on; grid minor

subplot(3,2,3); 
histogram(T2F_Picked_T0,'BinWidth',0.002,'DisplayStyle','stairs','LineWidth',2); hold on 
histogram(T2F_Picked_T1,'BinWidth',0.002,'DisplayStyle','stairs','LineWidth',2); hold on; xlabel('T_{2F} (s)','FontSize',18); ylabel('Count','FontSize',18); 
histogram(T2F_Picked_T2,'BinWidth',0.002,'DisplayStyle','stairs','LineWidth',2); 
histogram(T2F_Picked_T3,'BinWidth',0.002,'DisplayStyle','stairs','LineWidth',2); 
line([T2F_1 T2F_1],[0 260],'Color',cm(2,:),'LineStyle',':','LineWidth',4); line([T2F_2 T2F_2],[0 260],'Color',cm(3,:),'LineStyle',':','LineWidth',4);  line([T2F_3 T2F_3],[0 260],'Color',cm(4,:),'LineStyle',':','LineWidth',4); line([T2F_0 T2F_0],[0 260],'Color',cm(1,:),'LineStyle',':','LineWidth',4)
get(gca, 'XTick'); set(gca, 'FontSize', 18); get(gca, 'YTick'); set(gca, 'FontSize', 18); ylim([0 260]);
grid on; grid minor

subplot(3,2,4); 
histogram(T2S_Picked_T0,'BinWidth',0.005,'DisplayStyle','stairs','LineWidth',2); hold on 
histogram(T2S_Picked_T1,'BinWidth',0.005,'DisplayStyle','stairs','LineWidth',2); hold on; xlabel('T_{2S} (s)','FontSize',18); ylabel('Count','FontSize',18); 
histogram(T2S_Picked_T2,'BinWidth',0.005,'DisplayStyle','stairs','LineWidth',2); 
histogram(T2S_Picked_T3,'BinWidth',0.005,'DisplayStyle','stairs','LineWidth',2); 
line([T2S_1 T2S_1],[0 150],'Color',cm(2,:),'LineStyle',':','LineWidth',4); line([T2S_2 T2S_2],[0 150],'Color',cm(3,:),'LineStyle',':','LineWidth',4);  line([T2S_3 T2S_3],[0 150],'Color',cm(4,:),'LineStyle',':','LineWidth',4); line([T2S_0 T2S_0],[0 150],'Color',cm(1,:),'LineStyle',':','LineWidth',4)
get(gca, 'XTick'); set(gca, 'FontSize', 18); get(gca, 'YTick'); set(gca, 'FontSize', 18); ylim([0 150]);
grid on; grid minor

subplot(3,2,5); 
histogram(M0F_Picked_T0,'BinWidth',0.02,'DisplayStyle','stairs','LineWidth',2); hold on 
histogram(M0F_Picked_T1,'BinWidth',0.02,'DisplayStyle','stairs','LineWidth',2); hold on; xlabel('MWF','FontSize',18); ylabel('Count','FontSize',18); 
histogram(M0F_Picked_T2,'BinWidth',0.02,'DisplayStyle','stairs','LineWidth',2); 
histogram(M0F_Picked_T3,'BinWidth',0.02,'DisplayStyle','stairs','LineWidth',2); 
line([M0F_1 M0F_1],[0 200],'Color',cm(2,:),'LineStyle',':','LineWidth',4); line([M0F_2 M0F_2],[0 200],'Color',cm(3,:),'LineStyle',':','LineWidth',4);  line([M0F_3 M0F_3],[0 200],'Color',cm(4,:),'LineStyle',':','LineWidth',4); line([M0F_0 M0F_0],[0 200],'Color',cm(1,:),'LineStyle',':','LineWidth',4)
get(gca, 'XTick'); set(gca, 'FontSize', 18); get(gca, 'YTick'); set(gca, 'FontSize', 18);
grid on; grid minor

switch exchange
    case 'on'
        subplot(3,2,6);
        histogram(kFS_Picked_T0,'BinWidth',2,'DisplayStyle','stairs','LineWidth',2); hold on
        histogram(kFS_Picked_T1,'BinWidth',2,'DisplayStyle','stairs','LineWidth',2); hold on; xlabel('k_{FS} (s^{-1})','FontSize',18); ylabel('Count','FontSize',18);
        histogram(kFS_Picked_T2,'BinWidth',2,'DisplayStyle','stairs','LineWidth',2);
        histogram(kFS_Picked_T3,'BinWidth',2,'DisplayStyle','stairs','LineWidth',2);
        line([kFS_1 kFS_1],[0 200],'Color',cm(2,:),'LineStyle',':','LineWidth',4); line([kFS_2 kFS_2],[0 200],'Color',cm(3,:),'LineStyle',':','LineWidth',4);  line([kFS_3 kFS_3],[0 200],'Color',cm(4,:),'LineStyle',':','LineWidth',4); line([kFS_0 kFS_0],[0 200],'Color',cm(1,:),'LineStyle',':','LineWidth',4)
        get(gca, 'XTick'); set(gca, 'FontSize', 18); get(gca, 'YTick'); set(gca, 'FontSize', 18); ylim([0 200])
        grid on; grid minor
    case 'off'
end

%% Degenerate signal plots.

SPGR_Data_0 = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1S_0,'T1_F',T1F_0,'M0_F',M0F_0,'k_FS',kFS_0);
SSFP_Data_0_0 = SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1_0,'T1_S',T1S_0,'T2_S',T2S_0,'T1_F',T1F_0,'T2_F',T2F_0,'M0_F',M0F_0,'k_FS',kFS_0);
SSFP_Data_180_0 = SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2_0,'T1_S',T1S_0,'T2_S',T2S_0,'T1_F',T1F_0,'T2_F',T2F_0,'M0_F',M0F_0,'k_FS',kFS_0);
SSFP_Data_0 = [SSFP_Data_0_0 ; SSFP_Data_180_0];

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

figure(2); subplot(2,2,1)
plot(rad2deg(FA_SPGR), SPGR_Data_0/mean(SPGR_Data_0), 'mo','Linewidth',2,'MarkerSize',10); hold on
plot(rad2deg(FA_SSFP0), SSFP_Data_0_0/mean(SSFP_Data_0), 'ko','LineWidth',2,'MarkerSize',10)
plot(rad2deg(FA_SSFP180), SSFP_Data_180_0/mean(SSFP_Data_0), 'ro','LineWidth',2,'MarkerSize',10)
xlabel('FA (^{o})', 'FontSize', 18); ylabel('Normalised Signal (a.u.)', 'FontSize',18);
ll = legend({'SPGR','bSSFP0','bSSFP180'},'Orientation','horizontal','Position',[0.395810568652952,0.490026486011197,0.214207646283296,0.039999998966853]); ll.FontSize = 20; ll.AutoUpdate = 'off'; legend('boxoff'); 
get(gca, 'XTick'); set(gca, 'FontSize', 18); get(gca, 'YTick'); set(gca, 'FontSize', 18); tt = title('HB'); tt.FontSize = 24; grid on; grid minor;
subplot(2,2,2)
plot(rad2deg(FA_SPGR), SPGR_Data_1/mean(SPGR_Data_1), 'mo','Linewidth',2,'MarkerSize',10); hold on
plot(rad2deg(FA_SSFP0), SSFP_Data_0_1/mean(SSFP_Data_1), 'ko','LineWidth',2,'MarkerSize',10)
plot(rad2deg(FA_SSFP180), SSFP_Data_180_1/mean(SSFP_Data_1), 'ro','LineWidth',2,'MarkerSize',10)
xlabel('FA (^{o})', 'FontSize', 18); ylabel('Normalised Signal (a.u.)', 'FontSize',18);
get(gca, 'XTick'); set(gca, 'FontSize', 18); get(gca, 'YTick'); set(gca, 'FontSize', 18); tt = title('WML'); tt.FontSize = 24; grid on; grid minor;
subplot(2,2,3)
plot(rad2deg(FA_SPGR), SPGR_Data_2/mean(SPGR_Data_2), 'mo','Linewidth',2,'MarkerSize',10); hold on
plot(rad2deg(FA_SSFP0), SSFP_Data_0_2/mean(SSFP_Data_2), 'ko','LineWidth',2,'MarkerSize',10)
plot(rad2deg(FA_SSFP180), SSFP_Data_180_2/mean(SSFP_Data_2), 'ro','LineWidth',2,'MarkerSize',10)
xlabel('FA (^{o})', 'FontSize', 18); ylabel('Normalised Signal (a.u.)', 'FontSize',18);
get(gca, 'XTick'); set(gca, 'FontSize', 18); get(gca, 'YTick'); set(gca, 'FontSize', 18); tt = title('INT'); tt.FontSize = 24; grid on; grid minor;
subplot(2,2,4)
plot(rad2deg(FA_SPGR), SPGR_Data_3/mean(SPGR_Data_3), 'mo','Linewidth',2,'MarkerSize',10); hold on
plot(rad2deg(FA_SSFP0), SSFP_Data_0_3/mean(SSFP_Data_3), 'ko','LineWidth',2,'MarkerSize',10)
plot(rad2deg(FA_SSFP180), SSFP_Data_180_3/mean(SSFP_Data_3), 'ro','LineWidth',2,'MarkerSize',10)
xlabel('FA (^{o})', 'FontSize', 18); ylabel('Normalised Signal (a.u.)', 'FontSize',18);
get(gca, 'XTick'); set(gca, 'FontSize', 18); get(gca, 'YTick'); set(gca, 'FontSize', 18); tt = title('GML'); tt.FontSize = 24; grid on; grid minor;

SPGR_FM_0 = zeros(length(FA_SPGR),PlottingNo);
SSFP0_FM_0 = zeros(length(FA_SSFP0),PlottingNo);
SSFP180_FM_0 = zeros(length(FA_SSFP180),PlottingNo);
SSFP_FM_0 = zeros((length(FA_SSFP0)+length(FA_SSFP180)),PlottingNo);
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
    
    subplot(2,2,1)
    SPGR_FM_0(:,ii) = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1S_Picked_T0(ii),'T1_F',T1F_Picked_T0(ii),'M0_F',M0F_Picked_T0(ii),'k_FS',kFS_Picked_T0(ii));
    SSFP0_FM_0(:,ii) = SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1_0,'T1_S',T1S_Picked_T0(ii),'T2_S',T2S_Picked_T0(ii),'T1_F',T1F_Picked_T0(ii),'T2_F',T2F_Picked_T0(ii),'M0_F',M0F_Picked_T0(ii),'k_FS',kFS_Picked_T0(ii));
    SSFP180_FM_0(:,ii) = SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2_0,'T1_S',T1S_Picked_T0(ii),'T2_S',T2S_Picked_T0(ii),'T1_F',T1F_Picked_T0(ii),'T2_F',T2F_Picked_T0(ii),'M0_F',M0F_Picked_T0(ii),'k_FS',kFS_Picked_T0(ii));
    SSFP_FM_0(:,ii) = [SSFP0_FM_0(:,ii) ; SSFP180_FM_0(:,ii)];    
    subplot(2,2,2)
    SPGR_FM_1(:,ii) = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1S_Picked_T1(ii),'T1_F',T1F_Picked_T1(ii),'M0_F',M0F_Picked_T1(ii),'k_FS',kFS_Picked_T1(ii));
    SSFP0_FM_1(:,ii) = SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1_1,'T1_S',T1S_Picked_T1(ii),'T2_S',T2S_Picked_T1(ii),'T1_F',T1F_Picked_T1(ii),'T2_F',T2F_Picked_T1(ii),'M0_F',M0F_Picked_T1(ii),'k_FS',kFS_Picked_T1(ii));
    SSFP180_FM_1(:,ii) = SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2_1,'T1_S',T1S_Picked_T1(ii),'T2_S',T2S_Picked_T1(ii),'T1_F',T1F_Picked_T1(ii),'T2_F',T2F_Picked_T1(ii),'M0_F',M0F_Picked_T1(ii),'k_FS',kFS_Picked_T1(ii));
    SSFP_FM_1(:,ii) = [SSFP0_FM_1(:,ii) ; SSFP180_FM_1(:,ii)];
    subplot(2,2,3)
    SPGR_FM_2(:,ii) = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1S_Picked_T2(ii),'T1_F',T1F_Picked_T2(ii),'M0_F',M0F_Picked_T2(ii),'k_FS',kFS_Picked_T2(ii));
    SSFP0_FM_2(:,ii) = SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1_2,'T1_S',T1S_Picked_T2(ii),'T2_S',T2S_Picked_T2(ii),'T1_F',T1F_Picked_T2(ii),'T2_F',T2F_Picked_T2(ii),'M0_F',M0F_Picked_T2(ii),'k_FS',kFS_Picked_T2(ii));
    SSFP180_FM_2(:,ii) = SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2_2,'T1_S',T1S_Picked_T2(ii),'T2_S',T2S_Picked_T2(ii),'T1_F',T1F_Picked_T2(ii),'T2_F',T2F_Picked_T2(ii),'M0_F',M0F_Picked_T2(ii),'k_FS',kFS_Picked_T2(ii));
    SSFP_FM_2(:,ii) = [SSFP0_FM_2(:,ii) ; SSFP180_FM_2(:,ii)];
    subplot(2,2,4)
    SPGR_FM_3(:,ii) = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1S_Picked_T3(ii),'T1_F',T1F_Picked_T3(ii),'M0_F',M0F_Picked_T3(ii),'k_FS',kFS_Picked_T3(ii));
    SSFP0_FM_3(:,ii) = SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1_3,'T1_S',T1S_Picked_T3(ii),'T2_S',T2S_Picked_T3(ii),'T1_F',T1F_Picked_T3(ii),'T2_F',T2F_Picked_T3(ii),'M0_F',M0F_Picked_T3(ii),'k_FS',kFS_Picked_T3(ii));
    SSFP180_FM_3(:,ii) = SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2_3,'T1_S',T1S_Picked_T3(ii),'T2_S',T2S_Picked_T3(ii),'T1_F',T1F_Picked_T3(ii),'T2_F',T2F_Picked_T3(ii),'M0_F',M0F_Picked_T3(ii),'k_FS',kFS_Picked_T3(ii));
    SSFP_FM_3(:,ii) = [SSFP0_FM_3(:,ii) ; SSFP180_FM_3(:,ii)];
    
    SPGR_FM_0(:,ii) = SPGR_FM_0(:,ii)/mean(SPGR_FM_0(:,ii));
    SSFP0_FM_0(:,ii) = SSFP0_FM_0(:,ii)/mean(SSFP_FM_0(:,ii));
    SSFP180_FM_0(:,ii) = SSFP180_FM_0(:,ii)/mean(SSFP_FM_0(:,ii));
    
    SPGR_FM_1(:,ii) = SPGR_FM_1(:,ii)/mean(SPGR_FM_1(:,ii));
    SSFP0_FM_1(:,ii) = SSFP0_FM_1(:,ii)/mean(SSFP_FM_1(:,ii));
    SSFP180_FM_1(:,ii) = SSFP180_FM_1(:,ii)/mean(SSFP_FM_1(:,ii));
    
    SPGR_FM_2(:,ii) = SPGR_FM_2(:,ii)/mean(SPGR_FM_2(:,ii));
    SSFP0_FM_2(:,ii) = SSFP0_FM_2(:,ii)/mean(SSFP_FM_2(:,ii));
    SSFP180_FM_2(:,ii) = SSFP180_FM_2(:,ii)/mean(SSFP_FM_2(:,ii));
    
    SPGR_FM_3(:,ii) = SPGR_FM_3(:,ii)/mean(SPGR_FM_3(:,ii));
    SSFP0_FM_3(:,ii) = SSFP0_FM_3(:,ii)/mean(SSFP_FM_3(:,ii));
    SSFP180_FM_3(:,ii) = SSFP180_FM_3(:,ii)/mean(SSFP_FM_3(:,ii));
    
end

subplot(2,2,1)
plot(rad2deg(FA_SPGR), SPGR_FM_0, '-', 'Color', cm(1,:), 'LineWidth',0.1);
plot(rad2deg(FA_SSFP0), SSFP0_FM_0, '-', 'Color', cm(1,:), 'LineWidth',0.1);
plot(rad2deg(FA_SSFP180), SSFP180_FM_0, '-', 'Color', cm(1,:), 'LineWidth',0.1); hold on;
subplot(2,2,2)
plot(rad2deg(FA_SPGR), SPGR_FM_1, '-', 'Color', cm(2,:), 'LineWidth',0.1);
plot(rad2deg(FA_SSFP0), SSFP0_FM_1, '-', 'Color', cm(2,:), 'LineWidth',0.1);
plot(rad2deg(FA_SSFP180), SSFP180_FM_1, '-', 'Color', cm(2,:), 'LineWidth',0.1); hold on;
subplot(2,2,3)
plot(rad2deg(FA_SPGR), SPGR_FM_2, '-', 'Color', cm(3,:), 'LineWidth',0.1);
plot(rad2deg(FA_SSFP0), SSFP0_FM_2, '-', 'Color', cm(3,:), 'LineWidth',0.1);
plot(rad2deg(FA_SSFP180), SSFP180_FM_2, '-', 'Color', cm(3,:), 'LineWidth',0.1); hold on;
subplot(2,2,4)
plot(rad2deg(FA_SPGR), SPGR_FM_3, '-', 'Color', cm(4,:), 'LineWidth',0.1);
plot(rad2deg(FA_SSFP0), SSFP0_FM_3, '-', 'Color', cm(4,:), 'LineWidth',0.1);
plot(rad2deg(FA_SSFP180), SSFP180_FM_3, '-', 'Color', cm(4,:), 'LineWidth',0.1); hold on;

subplot(2,2,1)
hold on;
plot([29.6 33.5],[2.068 1.21],'--','LineWidth',1.5,'Color',[0.5 0.5 0.5])
plot([30.4 54.2],[2.1 1.21],'--','LineWidth',1.5,'Color',[0.5 0.5 0.5])
axes('Position',[.29 .65 .1 .1]); box on
hnd1 = plot(rad2deg(FA_SSFP180), SSFP180_FM_0, '-', 'Color', cm(1,:), 'LineWidth',0.1); hold on;
ylim([2.068 2.1]); xlim([29.6 30.4]);
get(gca, 'XTick'); set(gca, 'FontSize', 14); get(gca, 'YTick'); set(gca, 'FontSize', 14)