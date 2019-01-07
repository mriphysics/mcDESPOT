%% Format CSMT data.

% Concatenate Data.  
SPGR = cat(4, SPGRFA2_Image, SPGRFA4_Image, SPGRFA6_Image, SPGRFA8_Image, SPGRFA10_Image, SPGRFA12_Image, SPGRFA14_Image, SPGRFA16_Image, SPGRFA18_Image, SPGRFA20_Image);
SSFP = cat(4, SSFPFA2_Image, SSFPFA6_Image, SSFPFA14_Image, SSFPFA22_Image, SSFPFA30_Image, SSFPFA38_Image, SSFPFA46_Image, SSFPFA54_Image, SSFPFA62_Image, SSFPFA70_Image, SSFPFA14PC0_Image, SSFPFA22PC0_Image, SSFPFA30PC0_Image, SSFPFA38PC0_Image, SSFPFA46PC0_Image, SSFPFA54PC0_Image, SSFPFA62PC0_Image, SSFPFA70PC0_Image);

% SRC parameters.
Trials = 40000; Iterations = 30; N = 50;

% Acquisition parameters.
TR_SPGR = 6.5e-3; TR_SSFP = 6.5e-3; 
FA_SPGR_Original = [2 4 6 8 10 12 14 16 18 20];
FA_SSFP_Original = [2 6 14 22 30 38 46 54 62 70];
FA_SSFP0_Original = [14 22 30 38 46 54 62 70];

% Define region.
Slice = 64; x_min = 11; x_max = 121; y_min = 1; y_max = 138;
y_vector = (y_min:1:y_max)'; x_vector = (x_min:1:x_max)';
Coords = zeros(length(y_vector)*length(x_vector),3);

Step = 0;
for ii = 1:length(y_vector)
    for jj = 1:length(x_vector)
        Step = Step + 1;
        Coords(Step,:) = [Slice,y_vector(ii),x_vector(jj)];
    end
end

% Initialise loop variables.
SPGR_Data = zeros(length(Coords),length(FA_SPGR_Original)); SSFP_Data = zeros(length(Coords),(length(FA_SSFP_Original)+length(FA_SSFP0_Original)));
SPGR_Data_Norm = zeros(length(Coords),length(FA_SPGR_Original)); SSFP_Data_Norm = zeros(length(Coords),(length(FA_SSFP_Original)+length(FA_SSFP0_Original)));
B1_Data = zeros(length(Coords),1);

FA_SPGR_Rad = zeros(length(Coords),length(FA_SPGR_Original)); FA_SSFP180_Rad = zeros(length(Coords),length(FA_SSFP_Original)); FA_SSFP0_Rad = zeros(length(Coords),length(FA_SSFP0_Original));
FA_SPGR = zeros(length(Coords),length(FA_SPGR_Original)); FA_SSFP180 = zeros(length(Coords),length(FA_SSFP_Original)); FA_SSFP0 = zeros(length(Coords),length(FA_SSFP0_Original));
Data_Concatenated = zeros((length(FA_SPGR_Original)+length(FA_SSFP_Original)+length(FA_SSFP0_Original)),length(Coords));

for ll = 1:length(Coords)
    
    % Extract and concatenate ground-truth signals.
    SPGR_Data(ll,:) = abs(SPGR(Coords(ll,1),Coords(ll,2),Coords(ll,3),:));
    SSFP_Data(ll,:) = abs(SSFP(Coords(ll,1),Coords(ll,2),Coords(ll,3),:));
    
    SPGR_Data_Norm(ll,:) = SPGR_Data(ll,:)./mean(SPGR_Data(ll,:),2);
    SSFP_Data_Norm(ll,:) = SSFP_Data(ll,:)./mean(SSFP_Data(ll,:),2);
    
end
    
for ll = 1:length(Coords)
    
    Data_Concatenated(:,ll) = [SPGR_Data_Norm(ll,:), SSFP_Data_Norm(ll,:)].';
    
    % Correct FAs using B1-map.
    B1_Data(ll) = B1_Image(Coords(ll,1),Coords(ll,2),Coords(ll,3));
    FA_SPGR_Rad(ll,:) = B1_Data(ll) .* deg2rad(FA_SPGR_Original); FA_SPGR = rad2deg(FA_SPGR_Rad);
    FA_SSFP180_Rad(ll,:) = B1_Data(ll) .* deg2rad(FA_SSFP_Original); FA_SSFP180 = rad2deg(FA_SSFP180_Rad);
    FA_SSFP0_Rad(ll,:) = B1_Data(ll) .* deg2rad(FA_SSFP0_Original); FA_SSFP0 = rad2deg(FA_SSFP0_Rad);
    
end

% Only fit for brain pixels.
Indices_Fitted = find(SPGR_Data(:,10) > 1);
Fitted_Data = Data_Concatenated(:,Indices_Fitted); Fitted_Coords = Coords(Indices_Fitted,:);
Fitted_FA_SPGR = FA_SPGR(Indices_Fitted,:); Fitted_FA_SSFP180 = FA_SSFP180(Indices_Fitted,:); Fitted_FA_SSFP0 = FA_SSFP0(Indices_Fitted,:);

delete(gcp('nocreate')); c = parcluster('local'); c.NumWorkers = 8; parpool(c, c.NumWorkers);
[T1S_Bouhrara, T1F_Bouhrara, T2S_Bouhrara, T2F_Bouhrara, kFS_Bouhrara, M0F_Bouhrara, Delta_Bouhrara] = SRC_mcDESPOT_Bouhrara(Fitted_Coords, Trials, Iterations, N, Fitted_FA_SPGR, Fitted_FA_SSFP180, Fitted_FA_SSFP0, TR_SPGR, TR_SSFP, Fitted_Data);

MWF_Bouhrara = zeros(size(Coords,1),1);
Delta_Bouhrara_2 = zeros(size(Coords,1),1);
kFS_Bouhrara_2 = zeros(size(Coords,1),1);
T1F_Bouhrara_2 = zeros(size(Coords,1),1);
T1S_Bouhrara_2 = zeros(size(Coords,1),1);
T2F_Bouhrara_2 = zeros(size(Coords,1),1);
T2S_Bouhrara_2 = zeros(size(Coords,1),1);
for pp = 1:length(Indices_Fitted)
   
    MWF_Bouhrara(Indices_Fitted(pp),1) = M0F_Bouhrara(pp);
    Delta_Bouhrara_2(Indices_Fitted(pp),1) = Delta_Bouhrara(pp);
    kFS_Bouhrara_2(Indices_Fitted(pp),1) = kFS_Bouhrara(pp);
    T1F_Bouhrara_2(Indices_Fitted(pp),1) = T1F_Bouhrara(pp);
    T1S_Bouhrara_2(Indices_Fitted(pp),1) = T1S_Bouhrara(pp);
    T2F_Bouhrara_2(Indices_Fitted(pp),1) = T2F_Bouhrara(pp);
    T2S_Bouhrara_2(Indices_Fitted(pp),1) = T2S_Bouhrara(pp);

end

Map_Bouhrara = (vec2mat(MWF_Bouhrara,length(x_vector)));

%% Format non-CSMT data.

% Concatenate Data.  
SPGR = cat(4, SPGRFA2_Image, SPGRFA4_Image, SPGRFA6_Image, SPGRFA8_Image, SPGRFA10_Image, SPGRFA12_Image, SPGRFA14_Image, SPGRFA16_Image, SPGRFA18_Image, SPGRFA20_Image);
SSFP = cat(4, SSFPFA2_Image, SSFPFA6_Image, SSFPFA14_Image, SSFPFA22_Image, SSFPFA30_Image, SSFPFA38_Image, SSFPFA46_Image, SSFPFA54_Image, SSFPFA62_Image, SSFPFA70_Image, SSFPFA14PC0_Image, SSFPFA22PC0_Image, SSFPFA30PC0_Image, SSFPFA38PC0_Image, SSFPFA46PC0_Image, SSFPFA54PC0_Image, SSFPFA62PC0_Image, SSFPFA70PC0_Image);

% SRC parameters.
Trials = 40000; Iterations = 30; N = 50;

% Acquisition parameters.
TR_SPGR = 6.5e-3; TR_SSFP = 6.5e-3; 
FA_SPGR_Original = [2 4 6 8 10 12 14 16 18 20];
FA_SSFP_Original = [2 6 14 22 30 38 46 54 62 70];
FA_SSFP0_Original = [14 22 30 38 46 54 62 70];

% Define region.
Slice = 74; x_min = 11; x_max = 121; y_min = 1; y_max = 138;
y_vector = (y_min:1:y_max)'; x_vector = (x_min:1:x_max)';
Coords = zeros(length(y_vector)*length(x_vector),3);

Step = 0;
for ii = 1:length(y_vector)
    for jj = 1:length(x_vector)
        Step = Step + 1;
        Coords(Step,:) = [Slice,y_vector(ii),x_vector(jj)];
    end
end

% Initialise loop variables.
SPGR_Data = zeros(length(Coords),length(FA_SPGR_Original)); SSFP_Data = zeros(length(Coords),(length(FA_SSFP_Original)+length(FA_SSFP0_Original)));
SPGR_Data_Norm = zeros(length(Coords),length(FA_SPGR_Original)); SSFP_Data_Norm = zeros(length(Coords),(length(FA_SSFP_Original)+length(FA_SSFP0_Original)));
B1_Data = zeros(length(Coords),1);

FA_SPGR_Rad = zeros(length(Coords),length(FA_SPGR_Original)); FA_SSFP180_Rad = zeros(length(Coords),length(FA_SSFP_Original)); FA_SSFP0_Rad = zeros(length(Coords),length(FA_SSFP0_Original));
FA_SPGR = zeros(length(Coords),length(FA_SPGR_Original)); FA_SSFP180 = zeros(length(Coords),length(FA_SSFP_Original)); FA_SSFP0 = zeros(length(Coords),length(FA_SSFP0_Original));
Data_Concatenated = zeros((length(FA_SPGR_Original)+length(FA_SSFP_Original)+length(FA_SSFP0_Original)),length(Coords));

for ll = 1:length(Coords)
    
    % Extract and concatenate ground-truth signals.
    SPGR_Data(ll,:) = abs(SPGR(Coords(ll,1),Coords(ll,2),Coords(ll,3),:));
    SSFP_Data(ll,:) = abs(SSFP(Coords(ll,1),Coords(ll,2),Coords(ll,3),:));
    
    SPGR_Data_Norm(ll,:) = SPGR_Data(ll,:)./mean(SPGR_Data(ll,:),2);
    SSFP_Data_Norm(ll,:) = SSFP_Data(ll,:)./mean(SSFP_Data(ll,:),2);
    
end
    
for ll = 1:length(Coords)
    
    Data_Concatenated(:,ll) = [SPGR_Data_Norm(ll,:), SSFP_Data_Norm(ll,:)].';
    
    % Correct FAs using B1-map.
    B1_Data(ll) = B1_Image(Coords(ll,1),Coords(ll,2),Coords(ll,3));
    FA_SPGR_Rad(ll,:) = B1_Data(ll) .* deg2rad(FA_SPGR_Original); FA_SPGR = rad2deg(FA_SPGR_Rad);
    FA_SSFP180_Rad(ll,:) = B1_Data(ll) .* deg2rad(FA_SSFP_Original); FA_SSFP180 = rad2deg(FA_SSFP180_Rad);
    FA_SSFP0_Rad(ll,:) = B1_Data(ll) .* deg2rad(FA_SSFP0_Original); FA_SSFP0 = rad2deg(FA_SSFP0_Rad);
    
end

% Only fit for brain pixels.
Indices_Fitted = find(SPGR_Data(:,10) > 1);
Fitted_Data = Data_Concatenated(:,Indices_Fitted); Fitted_Coords = Coords(Indices_Fitted,:);
Fitted_FA_SPGR = FA_SPGR(Indices_Fitted,:); Fitted_FA_SSFP180 = FA_SSFP180(Indices_Fitted,:); Fitted_FA_SSFP0 = FA_SSFP0(Indices_Fitted,:);

delete(gcp('nocreate')); c = parcluster('local'); c.NumWorkers = 8; parpool(c, c.NumWorkers);
[NC_T1S_Bouhrara, NC_T1F_Bouhrara, NC_T2S_Bouhrara, NC_T2F_Bouhrara, NC_kFS_Bouhrara, NC_M0F_Bouhrara, NC_Delta_Bouhrara] = SRC_mcDESPOT_Bouhrara(Fitted_Coords, Trials, Iterations, N, Fitted_FA_SPGR, Fitted_FA_SSFP180, Fitted_FA_SSFP0, TR_SPGR, TR_SSFP, Fitted_Data);

NC_MWF_Bouhrara = zeros(size(NC_Coords,1),1);
NC_Delta_Bouhrara_2 = zeros(size(NC_Coords,1),1);
NC_kFS_Bouhrara_2 = zeros(size(NC_Coords,1),1);
NC_T1F_Bouhrara_2 = zeros(size(NC_Coords,1),1);
NC_T1S_Bouhrara_2 = zeros(size(NC_Coords,1),1);
NC_T2F_Bouhrara_2 = zeros(size(NC_Coords,1),1);
NC_T2S_Bouhrara_2 = zeros(size(NC_Coords,1),1);
for pp = 1:length(NC_Indices_Fitted)
   
    NC_MWF_Bouhrara(NC_Indices_Fitted(pp),1) = NC_M0F_Bouhrara(pp);
    NC_Delta_Bouhrara_2(NC_Indices_Fitted(pp),1) = NC_Delta_Bouhrara(pp);
    NC_kFS_Bouhrara_2(NC_Indices_Fitted(pp),1) = NC_kFS_Bouhrara(pp);
    NC_T1F_Bouhrara_2(NC_Indices_Fitted(pp),1) = NC_T1F_Bouhrara(pp);
    NC_T1S_Bouhrara_2(NC_Indices_Fitted(pp),1) = NC_T1S_Bouhrara(pp);
    NC_T2F_Bouhrara_2(NC_Indices_Fitted(pp),1) = NC_T2F_Bouhrara(pp);
    NC_T2S_Bouhrara_2(NC_Indices_Fitted(pp),1) = NC_T2S_Bouhrara(pp);

end

NC_Map_Bouhrara = (vec2mat(NC_MWF_Bouhrara,length(x_vector)));

%% Plots

% ROI = roipoly(Map_Bouhrara);
% MapCrop_Bouhrara = Map_Bouhrara; MapCrop_Bouhrara(~ROI) = 0;
% 
% ROI = roipoly(NC_Map_Bouhrara);
% NC_MapCrop_Bouhrara = Map_Bouhrara; MapCrop_Bouhrara(~ROI) = 0;

PixelNo = 8261; % 8171 = (68,74)

figure(1); subplot(2,5,1); imagesc(MapCrop_Bouhrara); cb1 = colorbar; axis off; axis square; colormap(magma); caxis([0 0.5]); cb1.FontSize = 12; hold on;

SPGR_GT = SPGR_Data(PixelNo,:)./(mean(SPGR_Data(PixelNo,:)));
SSFP_GT = SSFP_Data(PixelNo,:)./(mean(SSFP_Data(PixelNo,:)));
SSFP180_GT = SSFP_GT(1,1:10); SSFP0_GT = SSFP_GT(1,11:18);

[FM_SPGR_Bouhrara, FM_SSFP0_Bouhrara, FM_SSFP180_Bouhrara] = MEX_mcDESPOT_B0_DE_DiffFAs(T1F_Bouhrara_2(PixelNo), T1S_Bouhrara_2(PixelNo), T2F_Bouhrara_2(PixelNo), T2S_Bouhrara_2(PixelNo), kFS_Bouhrara_2(PixelNo), MWF_Bouhrara(PixelNo), TR_SPGR, TR_SSFP, Delta_Bouhrara_2(PixelNo), FA_SPGR(PixelNo,:), FA_SSFP0(PixelNo,:), FA_SSFP180(PixelNo,:));

FM_SSFP_Bouhrara = [FM_SSFP180_Bouhrara,FM_SSFP0_Bouhrara];
FM_SSFP_Bouhrara_Norm = FM_SSFP_Bouhrara./mean(FM_SSFP_Bouhrara);
FM_SPGR_Bouhrara_Norm = FM_SPGR_Bouhrara./mean(FM_SPGR_Bouhrara);
FM_SSFP180_Bouhrara_Norm = FM_SSFP_Bouhrara_Norm(1,1:10); FM_SSFP0_Bouhrara_Norm = FM_SSFP_Bouhrara_Norm(1,11:18);

figure(1); subplot(2,5,[2 3 4 5]);
plot(FA_SPGR(PixelNo,:),SPGR_GT,'bo','Linewidth',2,'MarkerSize',10); hold on
plot(FA_SSFP0(PixelNo,:),SSFP0_GT,'ko','Linewidth',2,'MarkerSize',10);
plot(FA_SSFP180(PixelNo,:),SSFP180_GT,'ro','Linewidth',2,'MarkerSize',10);
plot(FA_SPGR(PixelNo,:),FM_SPGR_Bouhrara_Norm,'Color','m','LineStyle','--','LineWidth',2)

ll = legend({'GT SPGR','GT bSSFP_{0}','GT bSSFP_{180}','Fitted Signal'}); ll.FontSize = 16; ll.AutoUpdate = 'off'; legend('boxoff');
title('CSMT','FontSize',16)

plot(FA_SSFP0(PixelNo,:),FM_SSFP0_Bouhrara_Norm,'Color','m','LineStyle','--','LineWidth',2); hold on
plot(FA_SSFP180(PixelNo,:),FM_SSFP180_Bouhrara_Norm,'Color','m','LineStyle','--','LineWidth',2); hold on

grid on; grid minor;
xlabel('FA (^o)','FontSize',14);ylabel('Normalised Signal (a.u.)','FontSize',14);
get(gca, 'XTick'); set(gca, 'FontSize', 14); get(gca, 'YTick'); set(gca, 'FontSize', 14);

figure(1); subplot(2,5,6); imagesc(NC_MapCrop_Bouhrara); cb1 = colorbar; axis off; axis square; colormap(magma); caxis([0 0.5]); cb1.FontSize = 12; hold on;

NC_SPGR_GT = NC_SPGR_Data(PixelNo,:)./(mean(NC_SPGR_Data(PixelNo,:)));
NC_SSFP_GT = NC_SSFP_Data(PixelNo,:)./(mean(NC_SSFP_Data(PixelNo,:)));
NC_SSFP180_GT = NC_SSFP_GT(1,1:10); NC_SSFP0_GT = SSFP_GT(1,11:18);

[NC_FM_SPGR_Bouhrara, NC_FM_SSFP0_Bouhrara, NC_FM_SSFP180_Bouhrara] = MEX_mcDESPOT_B0_DE_DiffFAs(NC_T1F_Bouhrara_2(PixelNo), NC_T1S_Bouhrara_2(PixelNo), NC_T2F_Bouhrara_2(PixelNo), NC_T2S_Bouhrara_2(PixelNo), NC_kFS_Bouhrara_2(PixelNo), NC_MWF_Bouhrara(PixelNo), NC_TR_SPGR, NC_TR_SSFP, NC_Delta_Bouhrara_2(PixelNo), NC_FA_SPGR(PixelNo,:), NC_FA_SSFP0(PixelNo,:), NC_FA_SSFP180(PixelNo,:));

NC_FM_SSFP_Bouhrara = [NC_FM_SSFP180_Bouhrara,NC_FM_SSFP0_Bouhrara];
NC_FM_SSFP_Bouhrara_Norm = NC_FM_SSFP_Bouhrara./mean(NC_FM_SSFP_Bouhrara);
NC_FM_SPGR_Bouhrara_Norm = NC_FM_SPGR_Bouhrara./mean(NC_FM_SPGR_Bouhrara);
NC_FM_SSFP180_Bouhrara_Norm = NC_FM_SSFP_Bouhrara_Norm(1,1:10); NC_FM_SSFP0_Bouhrara_Norm = NC_FM_SSFP_Bouhrara_Norm(1,11:18);

figure(1); subplot(2,5,[7 8 9 10]);
plot(NC_FA_SPGR(PixelNo,:),NC_SPGR_GT,'bo','Linewidth',2,'MarkerSize',10); hold on
plot(NC_FA_SSFP0(PixelNo,:),NC_SSFP0_GT,'ko','Linewidth',2,'MarkerSize',10);
plot(NC_FA_SSFP180(PixelNo,:),NC_SSFP180_GT,'ro','Linewidth',2,'MarkerSize',10);
plot(NC_FA_SPGR(PixelNo,:),NC_FM_SPGR_Bouhrara_Norm,'Color','m','LineStyle','--','LineWidth',2)

title('Non-CSMT','FontSize',16)

plot(NC_FA_SSFP0(PixelNo,:),NC_FM_SSFP0_Bouhrara_Norm,'Color','m','LineStyle','--','LineWidth',2); hold on
plot(NC_FA_SSFP180(PixelNo,:),NC_FM_SSFP180_Bouhrara_Norm,'Color','m','LineStyle','--','LineWidth',2); hold on

grid on; grid minor;
xlabel('FA (^o)','FontSize',14);ylabel('Normalised Signal (a.u.)','FontSize',14);
get(gca, 'XTick'); set(gca, 'FontSize', 14); get(gca, 'YTick'); set(gca, 'FontSize', 14);

%% Compute residuals.

CSMT_Data = [SPGR_Data./repmat(mean(SPGR_Data,2),[1 10]), SSFP_Data./repmat(mean(SSFP_Data,2),[1 18])];
nonCSMT_Data = [NC_SPGR_Data./repmat(mean(NC_SPGR_Data,2),[1 10]), NC_SSFP_Data./repmat(mean(NC_SSFP_Data,2),[1 18])];

FM_CSMT_Data = zeros(length(MWF_Bouhrara),28);
FM_nonCSMT_Data = zeros(length(NC_MWF_Bouhrara),28);
Residuals_CSMT = zeros(length(MWF_Bouhrara),1);
Residuals_nonCSMT = zeros(length(NC_MWF_Bouhrara),1);

for pp = 1:length(MWF_Bouhrara)
    
    [FM_SPGR_Bouhrara(pp,:), FM_SSFP0_Bouhrara(pp,:), FM_SSFP180_Bouhrara(pp,:)] = MEX_mcDESPOT_B0_DE_DiffFAs(T1F_Bouhrara_2(pp), T1S_Bouhrara_2(pp), T2F_Bouhrara_2(pp), T2S_Bouhrara_2(pp), kFS_Bouhrara_2(pp), MWF_Bouhrara(pp), TR_SPGR, TR_SSFP, Delta_Bouhrara_2(pp), FA_SPGR(pp,:), FA_SSFP0(pp,:), FA_SSFP180(pp,:));

    FM_SSFP_Bouhrara(pp,:) = [FM_SSFP180_Bouhrara(pp,:), FM_SSFP0_Bouhrara(pp,:)]; 
    
    FM_CSMT_Data(pp,:) = [FM_SPGR_Bouhrara(pp,:)./mean(FM_SPGR_Bouhrara(pp,:)), FM_SSFP_Bouhrara(pp,:)./mean(FM_SSFP_Bouhrara(pp,:))];

    Residuals_CSMT(pp) = norm(FM_CSMT_Data(pp,:)-CSMT_Data(pp,:));
    
end

for pp = 1:length(NC_MWF_Bouhrara)
    
    [NC_FM_SPGR_Bouhrara(pp,:), NC_FM_SSFP0_Bouhrara(pp,:), NC_FM_SSFP180_Bouhrara(pp,:)] = MEX_mcDESPOT_B0_DE_DiffFAs(NC_T1F_Bouhrara_2(pp), NC_T1S_Bouhrara_2(pp), NC_T2F_Bouhrara_2(pp), NC_T2S_Bouhrara_2(pp), NC_kFS_Bouhrara_2(pp), NC_MWF_Bouhrara(pp), NC_TR_SPGR, NC_TR_SSFP, NC_Delta_Bouhrara_2(pp), NC_FA_SPGR(pp,:), NC_FA_SSFP0(pp,:), NC_FA_SSFP180(pp,:));

    NC_FM_SSFP_Bouhrara(pp,:) = [NC_FM_SSFP180_Bouhrara(pp,:), NC_FM_SSFP0_Bouhrara(pp,:)]; 

    FM_nonCSMT_Data(pp,:) = [NC_FM_SPGR_Bouhrara(pp,:)./mean(NC_FM_SPGR_Bouhrara(pp,:)), NC_FM_SSFP_Bouhrara(pp,:)./mean(NC_FM_SSFP_Bouhrara(pp,:))];

    Residuals_nonCSMT(pp) = norm(FM_nonCSMT_Data(pp,:)-nonCSMT_Data(pp,:));
    
end

CSMT_ResidualMap = (vec2mat(Residuals_CSMT,length(x_vector)));
NC_ResidualMap = (vec2mat(Residuals_nonCSMT,length(NC_x_vector)));

% ROI = roipoly(CSMT_ResidualMap);
% CSMT_MapCrop_Residuals = CSMT_ResidualMap; CSMT_MapCrop_Residuals(~ROI) = 0;
% MapCrop_Bouhrara = Map_Bouhrara; MapCrop_Bouhrara(~ROI) = 0;
% 
% ROI = roipoly(NC_ResidualMap);
% NC_MapCrop_Residuals = NC_ResidualMap; NC_MapCrop_Residuals(~ROI) = 0;
% NC_MapCrop_Bouhrara = NC_Map_Bouhrara; NC_MapCrop_Bouhrara(~ROI) = 0;

figure(2); subplot(1,4,1); imagesc(MapCrop_Bouhrara); cb1 = colorbar; axis off; axis square; colormap(magma); caxis([0 0.5]); cb1.FontSize = 12; tt = title('CSMT MWF'); tt.FontSize = 16;
subplot(1,4,2); imagesc(CSMT_MapCrop_Residuals); cb1 = colorbar; axis off; axis square; colormap(magma); caxis([0 0.3]); cb1.FontSize = 12; tt = title('CSMT Residuals'); tt.FontSize = 16;
subplot(1,4,3); imagesc(NC_MapCrop_Bouhrara); cb1 = colorbar; axis off; axis square; colormap(magma); caxis([0 0.5]); cb1.FontSize = 12; tt = title('non-CSMT MWF'); tt.FontSize = 16;
subplot(1,4,4); imagesc(NC_MapCrop_Residuals); cb1 = colorbar; axis off; axis square; colormap(magma); caxis([0 0.3]); cb1.FontSize = 12; tt = title('non-CSMT Residuals'); tt.FontSize = 16;
