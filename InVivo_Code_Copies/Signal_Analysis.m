%%% Analyse in vivo signal fits for Figure 10. %%%

MWF_BouhraraRF = zeros(size(Coords,1),1);
MWF_DeoniRF = zeros(size(Coords,1),1);
MWF_WoodRF = zeros(size(Coords,1),1);
MWF_ZhangRF = zeros(size(Coords,1),1);
MWF_BouhraraLSQ = zeros(size(Coords,1),1);

kFS_BouhraraRF = zeros(size(Coords,1),1);
kFS_DeoniRF = zeros(size(Coords,1),1);
kFS_WoodRF = zeros(size(Coords,1),1);
kFS_ZhangRF = zeros(size(Coords,1),1);
kFS_BouhraraLSQ = zeros(size(Coords,1),1);

T1F_BouhraraRF = zeros(size(Coords,1),1);
T1F_DeoniRF = zeros(size(Coords,1),1);
T1F_WoodRF = zeros(size(Coords,1),1);
T1F_ZhangRF = zeros(size(Coords,1),1);
T1F_BouhraraLSQ = zeros(size(Coords,1),1);

T1S_BouhraraRF = zeros(size(Coords,1),1);
T1S_DeoniRF = zeros(size(Coords,1),1);
T1S_WoodRF = zeros(size(Coords,1),1);
T1S_ZhangRF = zeros(size(Coords,1),1);
T1S_BouhraraLSQ = zeros(size(Coords,1),1);

T2F_BouhraraRF = zeros(size(Coords,1),1);
T2F_DeoniRF = zeros(size(Coords,1),1);
T2F_WoodRF = zeros(size(Coords,1),1);
T2F_ZhangRF = zeros(size(Coords,1),1);
T2F_BouhraraLSQ = zeros(size(Coords,1),1);

T2S_BouhraraRF = zeros(size(Coords,1),1);
T2S_DeoniRF = zeros(size(Coords,1),1);
T2S_WoodRF = zeros(size(Coords,1),1);
T2S_ZhangRF = zeros(size(Coords,1),1);
T2S_BouhraraLSQ = zeros(size(Coords,1),1);

Delta_BouhraraRF = zeros(size(Coords,1),1);
Delta_DeoniRF = zeros(size(Coords,1),1);
Delta_WoodRF = zeros(size(Coords,1),1);
Delta_ZhangRF = zeros(size(Coords,1),1);
Delta_BouhraraLSQ = zeros(size(Coords,1),1);

for pp = 1:length(Indices_Fitted)
   
    MWF_BouhraraRF(Indices_Fitted(pp),1) = M0F_Bouhrara(pp);
    MWF_DeoniRF(Indices_Fitted(pp),1) = M0F_Deoni(pp);
    MWF_WoodRF(Indices_Fitted(pp),1) = M0F_Wood(pp);
    MWF_ZhangRF(Indices_Fitted(pp),1) = M0F_Zhang(pp);
    MWF_BouhraraLSQ(Indices_Fitted(pp),1) = Params_Est(pp,6);
    
    T1F_BouhraraRF(Indices_Fitted(pp),1) = T1F_Bouhrara(pp);
    T1F_DeoniRF(Indices_Fitted(pp),1) = T1F_Deoni(pp);
    T1F_WoodRF(Indices_Fitted(pp),1) = T1F_Wood(pp);
    T1F_ZhangRF(Indices_Fitted(pp),1) = T1F_Zhang(pp);
    T1F_BouhraraLSQ(Indices_Fitted(pp),1) = Params_Est(pp,1);
    
    T1S_BouhraraRF(Indices_Fitted(pp),1) = T1S_Bouhrara(pp);
    T1S_DeoniRF(Indices_Fitted(pp),1) = T1S_Deoni(pp);
    T1S_WoodRF(Indices_Fitted(pp),1) = T1S_Wood(pp);
    T1S_ZhangRF(Indices_Fitted(pp),1) = T1S_Zhang(pp);
    T1S_BouhraraLSQ(Indices_Fitted(pp),1) = Params_Est(pp,2);

    T2F_BouhraraRF(Indices_Fitted(pp),1) = T2F_Bouhrara(pp);
    T2F_DeoniRF(Indices_Fitted(pp),1) = T2F_Deoni(pp);
    T2F_WoodRF(Indices_Fitted(pp),1) = T2F_Wood(pp);
    T2F_ZhangRF(Indices_Fitted(pp),1) = T2F_Zhang(pp);
    T2F_BouhraraLSQ(Indices_Fitted(pp),1) = Params_Est(pp,3);  
    
    T2S_BouhraraRF(Indices_Fitted(pp),1) = T2S_Bouhrara(pp);
    T2S_DeoniRF(Indices_Fitted(pp),1) = T2S_Deoni(pp);
    T2S_WoodRF(Indices_Fitted(pp),1) = T2S_Wood(pp);
    T2S_ZhangRF(Indices_Fitted(pp),1) = T2S_Zhang(pp);
    T2S_BouhraraLSQ(Indices_Fitted(pp),1) = Params_Est(pp,4);     

    kFS_BouhraraRF(Indices_Fitted(pp),1) = kFS_Bouhrara(pp);
    kFS_DeoniRF(Indices_Fitted(pp),1) = kFS_Deoni(pp);
    kFS_WoodRF(Indices_Fitted(pp),1) = kFS_Wood(pp);
    kFS_ZhangRF(Indices_Fitted(pp),1) = kFS_Zhang(pp);
    kFS_BouhraraLSQ(Indices_Fitted(pp),1) = Params_Est(pp,5);  
    
    Delta_BouhraraRF(Indices_Fitted(pp),1) = Delta_Bouhrara(pp);
    Delta_DeoniRF(Indices_Fitted(pp),1) = Delta_Deoni(pp);
    Delta_WoodRF(Indices_Fitted(pp),1) = Delta_Wood(pp);
    Delta_ZhangRF(Indices_Fitted(pp),1) = Delta_Zhang(pp);
    Delta_BouhraraLSQ(Indices_Fitted(pp),1) = Params_Est(pp,7);      

end

%% Forward-model signals for chosen pixel and each bound set.

PixelNo = 9058; 

SPGR_GT = SPGR_Data(PixelNo,:)./(mean(SPGR_Data(PixelNo,:)));
SSFP_GT = SSFP_Data(PixelNo,:)./(mean(SSFP_Data(PixelNo,:)));

SSFP180_GT = SSFP_GT(1,1:10); SSFP0_GT = SSFP_GT(1,11:18);

[FM_SPGR_Bouhrara, FM_SSFP0_Bouhrara, FM_SSFP180_Bouhrara] = MEX_mcDESPOT_B0_DE_DiffFAs(T1F_BouhraraRF(PixelNo), T1S_BouhraraRF(PixelNo), T2F_BouhraraRF(PixelNo), T2S_BouhraraRF(PixelNo), kFS_BouhraraRF(PixelNo), MWF_BouhraraRF(PixelNo), TR_SPGR, TR_SSFP, Delta_BouhraraRF(PixelNo), FA_SPGR(PixelNo,:), FA_SSFP0(PixelNo,:), FA_SSFP180(PixelNo,:));
[FM_SPGR_Deoni, FM_SSFP0_Deoni, FM_SSFP180_Deoni] = MEX_mcDESPOT_B0_DE_DiffFAs(T1F_DeoniRF(PixelNo), T1S_DeoniRF(PixelNo), T2F_DeoniRF(PixelNo), T2S_DeoniRF(PixelNo), kFS_DeoniRF(PixelNo), MWF_DeoniRF(PixelNo), TR_SPGR, TR_SSFP, Delta_DeoniRF(PixelNo), FA_SPGR(PixelNo,:), FA_SSFP0(PixelNo,:), FA_SSFP180(PixelNo,:));
[FM_SPGR_Wood, FM_SSFP0_Wood, FM_SSFP180_Wood] = MEX_mcDESPOT_B0_DE_DiffFAs(T1F_WoodRF(PixelNo), T1S_WoodRF(PixelNo), T2F_WoodRF(PixelNo), T2S_WoodRF(PixelNo), kFS_WoodRF(PixelNo), MWF_WoodRF(PixelNo), TR_SPGR, TR_SSFP, Delta_WoodRF(PixelNo), FA_SPGR(PixelNo,:), FA_SSFP0(PixelNo,:), FA_SSFP180(PixelNo,:));
[FM_SPGR_Zhang, FM_SSFP0_Zhang, FM_SSFP180_Zhang] = MEX_mcDESPOT_B0_DE_DiffFAs(T1F_ZhangRF(PixelNo), T1S_ZhangRF(PixelNo), T2F_ZhangRF(PixelNo), T2S_ZhangRF(PixelNo), kFS_ZhangRF(PixelNo), MWF_ZhangRF(PixelNo), TR_SPGR, TR_SSFP, Delta_ZhangRF(PixelNo), FA_SPGR(PixelNo,:), FA_SSFP0(PixelNo,:), FA_SSFP180(PixelNo,:));
[FM_SPGR_BouhraraLSQ, FM_SSFP0_BouhraraLSQ, FM_SSFP180_BouhraraLSQ] = MEX_mcDESPOT_B0_DE_DiffFAs(T1F_BouhraraLSQ(PixelNo), T1S_BouhraraLSQ(PixelNo), T2F_BouhraraLSQ(PixelNo), T2S_BouhraraLSQ(PixelNo), kFS_BouhraraLSQ(PixelNo), MWF_BouhraraLSQ(PixelNo), TR_SPGR, TR_SSFP, Delta_BouhraraLSQ(PixelNo), FA_SPGR(PixelNo,:), FA_SSFP0(PixelNo,:), FA_SSFP180(PixelNo,:));

FM_SSFP_Bouhrara = [FM_SSFP180_Bouhrara,FM_SSFP0_Bouhrara];
FM_SSFP_Deoni = [FM_SSFP180_Deoni,FM_SSFP0_Deoni];
FM_SSFP_Wood = [FM_SSFP180_Wood,FM_SSFP0_Wood];
FM_SSFP_Zhang = [FM_SSFP180_Zhang,FM_SSFP0_Zhang];
FM_SSFP_BouhraraLSQ = [FM_SSFP180_BouhraraLSQ,FM_SSFP0_BouhraraLSQ];

FM_SSFP_Bouhrara_Norm = FM_SSFP_Bouhrara./mean(FM_SSFP_Bouhrara);
FM_SSFP_Deoni_Norm = FM_SSFP_Deoni./mean(FM_SSFP_Deoni);
FM_SSFP_Wood_Norm = FM_SSFP_Wood./mean(FM_SSFP_Wood);
FM_SSFP_Zhang_Norm = FM_SSFP_Zhang./mean(FM_SSFP_Zhang);
FM_SSFP_Bouhrara_NormLSQ = FM_SSFP_BouhraraLSQ./mean(FM_SSFP_BouhraraLSQ);

FM_SPGR_Bouhrara_Norm = FM_SPGR_Bouhrara./mean(FM_SPGR_Bouhrara);
FM_SPGR_Deoni_Norm = FM_SPGR_Deoni./mean(FM_SPGR_Deoni);
FM_SPGR_Wood_Norm = FM_SPGR_Wood./mean(FM_SPGR_Wood);
FM_SPGR_Zhang_Norm = FM_SPGR_Zhang./mean(FM_SPGR_Zhang);
FM_SPGR_Bouhrara_NormLSQ = FM_SPGR_BouhraraLSQ./mean(FM_SPGR_BouhraraLSQ);

FM_SSFP180_Bouhrara_Norm = FM_SSFP_Bouhrara_Norm(1,1:10); FM_SSFP0_Bouhrara_Norm = FM_SSFP_Bouhrara_Norm(1,11:18);
FM_SSFP180_Deoni_Norm = FM_SSFP_Deoni_Norm(1,1:10); FM_SSFP0_Deoni_Norm = FM_SSFP_Deoni_Norm(1,11:18);
FM_SSFP180_Wood_Norm = FM_SSFP_Wood_Norm(1,1:10); FM_SSFP0_Wood_Norm = FM_SSFP_Wood_Norm(1,11:18);
FM_SSFP180_Zhang_Norm = FM_SSFP_Zhang_Norm(1,1:10); FM_SSFP0_Zhang_Norm = FM_SSFP_Zhang_Norm(1,11:18);
FM_SSFP180_Bouhrara_NormLSQ = FM_SSFP_Bouhrara_NormLSQ(1,1:10); FM_SSFP0_Bouhrara_NormLSQ = FM_SSFP_Bouhrara_NormLSQ(1,11:18);

cm = colormap(lines(5));

figure(2); subplot(2,4,[7 8]);
plot(FA_SPGR(PixelNo,:),SPGR_GT,'bo','Linewidth',2,'MarkerSize',10); hold on
plot(FA_SSFP0(PixelNo,:),SSFP0_GT,'ko','Linewidth',2,'MarkerSize',10); 
plot(FA_SSFP180(PixelNo,:),SSFP180_GT,'ro','Linewidth',2,'MarkerSize',10);
plot(FA_SPGR(PixelNo,:),FM_SPGR_Bouhrara_Norm,'Color',cm(1,:),'LineStyle','--','LineWidth',2)
plot(FA_SPGR(PixelNo,:),FM_SPGR_Deoni_Norm,'Color',cm(2,:),'LineStyle','--','LineWidth',2)
plot(FA_SPGR(PixelNo,:),FM_SPGR_Wood_Norm,'Color',cm(3,:),'LineStyle','--','LineWidth',2)
plot(FA_SPGR(PixelNo,:),FM_SPGR_Zhang_Norm,'Color',cm(4,:),'LineStyle','--','LineWidth',2)
plot(FA_SPGR(PixelNo,:),FM_SPGR_Bouhrara_NormLSQ,'Color',cm(5,:),'LineStyle','--','LineWidth',2)

ll = legend({'GT SPGR','GT bSSFP_{0}','GT bSSFP_{180}','B1 Fitted Signal','B2 Fitted Signal','B3 Fitted Signal','B4 Fitted Signal','lsqnonlin Fitted Signal'},'Position',[0.543898005183711,0.442834147048852,0.376502723413739,0.137789900598877]);
ll.FontSize = 18;
ll.AutoUpdate = 'off'; legend('boxoff'); ll.NumColumns = 3;

plot(FA_SSFP0(PixelNo,:),FM_SSFP0_Bouhrara_Norm,'Color',cm(1,:),'LineStyle','--','LineWidth',2); hold on
plot(FA_SSFP0(PixelNo,:),FM_SSFP0_Deoni_Norm,'Color',cm(2,:),'LineStyle','--','LineWidth',2)
plot(FA_SSFP0(PixelNo,:),FM_SSFP0_Wood_Norm,'Color',cm(3,:),'LineStyle','--','LineWidth',2)
plot(FA_SSFP0(PixelNo,:),FM_SSFP0_Zhang_Norm,'Color',cm(4,:),'LineStyle','--','LineWidth',2)
plot(FA_SSFP0(PixelNo,:),FM_SSFP0_Bouhrara_NormLSQ,'Color',cm(5,:),'LineStyle','--','LineWidth',2);

plot(FA_SSFP180(PixelNo,:),FM_SSFP180_Bouhrara_Norm,'Color',cm(1,:),'LineStyle','--','LineWidth',2); hold on
plot(FA_SSFP180(PixelNo,:),FM_SSFP180_Deoni_Norm,'Color',cm(2,:),'LineStyle','--','LineWidth',2)
plot(FA_SSFP180(PixelNo,:),FM_SSFP180_Wood_Norm,'Color',cm(3,:),'LineStyle','--','LineWidth',2)
plot(FA_SSFP180(PixelNo,:),FM_SSFP180_Zhang_Norm,'Color',cm(4,:),'LineStyle','--','LineWidth',2)
plot(FA_SSFP180(PixelNo,:),FM_SSFP180_Bouhrara_NormLSQ,'Color',cm(5,:),'LineStyle','--','LineWidth',2);

grid on; grid minor;
xlabel('FA (^o)','FontSize',16);ylabel('Normalised Signal (a.u.)','FontSize',16);
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
text(0.02,0.88,'(g)','Units','Normalized','VerticalAlignment','Bottom','FontSize',16)

T1F_Vals = [T1F_BouhraraRF(PixelNo);T1F_DeoniRF(PixelNo);T1F_WoodRF(PixelNo);T1F_ZhangRF(PixelNo);T1F_BouhraraLSQ(PixelNo)];
T1S_Vals = [T1S_BouhraraRF(PixelNo);T1S_DeoniRF(PixelNo);T1S_WoodRF(PixelNo);T1S_ZhangRF(PixelNo);T1S_BouhraraLSQ(PixelNo)];
T2F_Vals = [T2F_BouhraraRF(PixelNo);T2F_DeoniRF(PixelNo);T2F_WoodRF(PixelNo);T2F_ZhangRF(PixelNo);T2F_BouhraraLSQ(PixelNo)];
T2S_Vals = [T2S_BouhraraRF(PixelNo);T2S_DeoniRF(PixelNo);T2S_WoodRF(PixelNo);T2S_ZhangRF(PixelNo);T2S_BouhraraLSQ(PixelNo)];
MWF_Vals = [MWF_BouhraraRF(PixelNo);MWF_DeoniRF(PixelNo);MWF_WoodRF(PixelNo);MWF_ZhangRF(PixelNo);MWF_BouhraraLSQ(PixelNo)];
kFS_Vals = [kFS_BouhraraRF(PixelNo);kFS_DeoniRF(PixelNo);kFS_WoodRF(PixelNo);kFS_ZhangRF(PixelNo);kFS_BouhraraLSQ(PixelNo)];
