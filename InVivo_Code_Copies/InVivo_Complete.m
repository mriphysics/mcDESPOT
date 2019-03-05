%%% mcDESPOT in vivo fitting script with data processed in "InVivo_Processing.m". %%%

%% Format data. !!! WITHOUT LOW FA bSSFP0 DATA. !!!

exchange = 'on';

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
Slice = 64; x_min = 11; x_max = 121; y_min = 1; y_max = 138; % 74 for non-CSMT.
y_vector = (y_min:1:y_max)'; x_vector = (x_min:1:x_max)';
Coords = zeros(length(y_vector)*length(x_vector),3);

figure(1); imagesc(abs(squeeze(B1_Image(Slice,:,:)))); axis square;

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

%% Stochastic region contraction.

delete(gcp('nocreate')); c = parcluster('local'); c.NumWorkers = 8; parpool(c, c.NumWorkers);

switch exchange
    case 'on'
        [T1S_Bouhrara, T1F_Bouhrara, T2S_Bouhrara, T2F_Bouhrara, kFS_Bouhrara, M0F_Bouhrara, Delta_Bouhrara] = SRC_mcDESPOT_Bouhrara(Fitted_Coords, Trials, Iterations, N, Fitted_FA_SPGR, Fitted_FA_SSFP180, Fitted_FA_SSFP0, TR_SPGR, TR_SSFP, Fitted_Data);
        
        [T1S_Deoni, T1F_Deoni, T2S_Deoni, T2F_Deoni, kFS_Deoni, M0F_Deoni, Delta_Deoni] = SRC_mcDESPOT_Deoni(Fitted_Coords, Trials, Iterations, N, Fitted_FA_SPGR, Fitted_FA_SSFP180, Fitted_FA_SSFP0, TR_SPGR, TR_SSFP, Fitted_Data);
        
        [T1S_Wood, T1F_Wood, T2S_Wood, T2F_Wood, kFS_Wood, M0F_Wood, Delta_Wood] = SRC_mcDESPOT_Wood(Fitted_Coords, Trials, Iterations, N, Fitted_FA_SPGR, Fitted_FA_SSFP180, Fitted_FA_SSFP0, TR_SPGR, TR_SSFP, Fitted_Data);
        
        [T1S_Zhang, T1F_Zhang, T2S_Zhang, T2F_Zhang, kFS_Zhang, M0F_Zhang, Delta_Zhang] = SRC_mcDESPOT_Zhang(Fitted_Coords, Trials, Iterations, N, Fitted_FA_SPGR, Fitted_FA_SSFP180, Fitted_FA_SSFP0, TR_SPGR, TR_SSFP, Fitted_Data);
    
    case 'off'
        [T1S_Bouhrara_NE, T1F_Bouhrara_NE, T2S_Bouhrara_NE, T2F_Bouhrara_NE, M0F_Bouhrara_NE, Delta_Bouhrara_NE] = SRC_mcDESPOT_Bouhrara_NoEx(Fitted_Coords, Trials, Iterations, N, Fitted_FA_SPGR, Fitted_FA_SSFP180, Fitted_FA_SSFP0, TR_SPGR, TR_SSFP, Fitted_Data);
        
        [T1S_Deoni_NE, T1F_Deoni_NE, T2S_Deoni_NE, T2F_Deoni_NE, M0F_Deoni_NE, Delta_Deoni_NE] = SRC_mcDESPOT_Deoni_NoEx(Fitted_Coords, Trials, Iterations, N, Fitted_FA_SPGR, Fitted_FA_SSFP180, Fitted_FA_SSFP0, TR_SPGR, TR_SSFP, Fitted_Data);
        
        [T1S_Wood_NE, T1F_Wood_NE, T2S_Wood_NE, T2F_Wood_NE, M0F_Wood_NE, Delta_Wood_NE] = SRC_mcDESPOT_Wood_NoEx(Fitted_Coords, Trials, Iterations, N, Fitted_FA_SPGR, Fitted_FA_SSFP180, Fitted_FA_SSFP0, TR_SPGR, TR_SSFP, Fitted_Data);
        
        [T1S_Zhang_NE, T1F_Zhang_NE, T2S_Zhang_NE, T2F_Zhang_NE, M0F_Zhang_NE, Delta_Zhang_NE] = SRC_mcDESPOT_Zhang_NoEx(Fitted_Coords, Trials, Iterations, N, Fitted_FA_SPGR, Fitted_FA_SSFP180, Fitted_FA_SSFP0, TR_SPGR, TR_SSFP, Fitted_Data);
end

%% Plot MWF maps and WM histograms for Figure 10. Modify for SF5.

MWF_Bouhrara = zeros(size(Coords,1),1);
MWF_Deoni = zeros(size(Coords,1),1);
MWF_Wood = zeros(size(Coords,1),1);
MWF_Zhang = zeros(size(Coords,1),1);

for pp = 1:length(Indices_Fitted)
   
    MWF_Bouhrara(Indices_Fitted(pp),1) = M0F_Bouhrara_NE(pp);
    MWF_Deoni(Indices_Fitted(pp),1) = M0F_Deoni_NE(pp);
    MWF_Wood(Indices_Fitted(pp),1) = M0F_Wood_NE(pp);
    MWF_Zhang(Indices_Fitted(pp),1) = M0F_Zhang_NE(pp);

end

Map_Bouhrara = flipud(vec2mat(MWF_Bouhrara,length(x_vector)));
Map_Deoni = flipud(vec2mat(MWF_Deoni,length(x_vector)));
Map_Wood = flipud(vec2mat(MWF_Wood,length(x_vector)));
Map_Zhang = flipud(vec2mat(MWF_Zhang,length(x_vector)));

ROI = roipoly(Map_Bouhrara); 

MapCrop_Bouhrara = Map_Bouhrara; MapCrop_Bouhrara(~ROI) = 0;
figure(2); subplot(2,4,1); 
imagesc(MapCrop_Bouhrara); hold on; cb1 = colorbar; axis off; colormap(magma); tt = title('B1'); tt.FontSize = 18; pbaspect([1.11 1.38 1]); caxis([0 0.35]); cb1.FontSize = 14; plot(67,57,'.','MarkerSize',20,'Color','g')
text(-20,130,'(a)','FontSize',16)

MapCrop_Deoni = Map_Deoni; MapCrop_Deoni(~ROI) = 0;
figure(2); subplot(2,4,2); imagesc(MapCrop_Deoni); hold on; cb2 = colorbar; axis off; colormap(magma); tt = title('B2'); tt.FontSize = 18; pbaspect([1.11 1.38 1]); caxis([0 0.35]); cb2.FontSize = 14; plot(67,57,'.','MarkerSize',20,'Color','g') % Was: 73,65
text(-20,130,'(b)','FontSize',16)

MapCrop_Wood = Map_Wood; MapCrop_Wood(~ROI) = 0;
figure(2); subplot(2,4,3); imagesc(MapCrop_Wood); hold on; cb3 = colorbar; axis off; colormap(magma); tt = title('B3'); tt.FontSize = 18; pbaspect([1.11 1.38 1]); caxis([0 0.35]); cb3.FontSize = 14; plot(67,57,'.','MarkerSize',20,'Color','g')
text(-20,130,'(c)','FontSize',16)

MapCrop_Zhang = Map_Zhang; MapCrop_Zhang(~ROI) = 0;
figure(2); subplot(2,4,4); imagesc(MapCrop_Zhang); hold on; cb4 = colorbar; axis off; colormap(magma); tt = title('B4'); tt.FontSize = 18; pbaspect([1.11 1.38 1]); caxis([0 0.35]); cb4.FontSize = 14; plot(67,57,'.','MarkerSize',20,'Color','g')
text(-20,130,'(d)','FontSize',16)

MapCrop_BouhraraLSQ = Map_BouhraraLSQ; MapCrop_BouhraraLSQ(~ROI) = 0; % From LM_Test.m.
figure(2); subplot(2,5,5); imagesc(MapCrop_BouhraraLSQ); hold on; cb1 = colorbar; axis off; colormap(magma); tt = title('lsqnonlin'); tt.FontSize = 18; pbaspect([1.11 1.38 1]); caxis([0 0.35]); cb1.FontSize = 14; plot(67,57,'.','MarkerSize',20,'Color','g')
text(-20,130,'(e)','FontSize',16)

WM_Data = zeros(length(Coords),1);
for ll = 1:length(Coords)
    WM_Data(ll,:) = WM_Mask_Image(Coords(ll,1),Coords(ll,2),Coords(ll,3),:);
end

WM_Data_Fitted = WM_Data(Indices_Fitted,:);

WM_Reformat = zeros(size(Coords,1),1);
for pp = 1:length(Indices_Fitted)
    WM_Reformat(Indices_Fitted(pp),1) = WM_Data_Fitted(pp);
end

WM_Map = flipud(vec2mat(WM_Reformat,length(x_vector)));

WMMap_Crop = WM_Map; WMMap_Crop(~ROI) = 0;

WM_Pixels = find(WMMap_Crop == 1);

figure(2); subplot(2,4,[5 6 7 8]);
[N1,E1] = histcounts(MapCrop_Bouhrara(WM_Pixels)); 
plot(E1(2:end),N1,'o--','Color',[0	0.447000000000000	0.741000000000000],'LineWidth',2); hold on
[N2,E2] = histcounts(MapCrop_Deoni(WM_Pixels));
plot(E2(2:end),N2,'o--','Color',[0.850000000000000	0.325000000000000	0.0980000000000000],'LineWidth',2)
[N3,E3] = histcounts(MapCrop_Wood(WM_Pixels)); 
plot(E3(2:end),N3,'o--','Color',[0.929000000000000	0.694000000000000	0.125000000000000],'LineWidth',2)
[N4,E4] = histcounts(MapCrop_Zhang(WM_Pixels)); 
plot(E4(2:end),N4,'o--','Color',[0.494000000000000	0.184000000000000	0.556000000000000],'LineWidth',2)
[N5,E5] = histcounts(MapCrop_BouhraraLSQ(WM_Pixels)); 
plot(E5(2:end),N5,'o--','Color',[0.466000000000000	0.674000000000000	0.188000000000000],'LineWidth',2)
xlabel('MWF','FontSize',18); ylabel('Count','FontSize',18)
ll = legend('B1','B2','B3','B4','lsqnonlin'); ll.FontSize = 18; legend('boxoff');
get(gca, 'XTick'); set(gca, 'FontSize', 14); get(gca, 'YTick'); set(gca, 'FontSize', 16);
grid on; grid minor;
text(0.02,0.88,'(f)','Units','Normalized','VerticalAlignment','Bottom','FontSize',16)

%% Re-format estimates and create maps for all parameters. Used for SF7.

T1F_MapVector = zeros(size(Coords,1),1); T1S_MapVector = zeros(size(Coords,1),1);
T2F_MapVector = zeros(size(Coords,1),1); T2S_MapVector = zeros(size(Coords,1),1);
kFS_MapVector = zeros(size(Coords,1),1); M0F_MapVector = zeros(size(Coords,1),1); 

for pp = 1:length(Indices_Fitted)
   
    T1F_MapVector(Indices_Fitted(pp),1) = T1F_Bouhrara(pp); T1S_MapVector(Indices_Fitted(pp),1) = T1S_Bouhrara(pp);
    T2F_MapVector(Indices_Fitted(pp),1) = T2F_Bouhrara(pp); T2S_MapVector(Indices_Fitted(pp),1) = T2S_Bouhrara(pp);
    kFS_MapVector(Indices_Fitted(pp),1) = kFS_Bouhrara(pp); M0F_MapVector(Indices_Fitted(pp),1) = M0F_Bouhrara(pp); 

end

M0F_Map = vec2mat(M0F_MapVector,length(x_vector));

Select ROI to discard erroneous regions at brain edges.
ROI = roipoly(M0F_Map); 

T1F_Map = vec2mat(T1F_MapVector,length(x_vector));
T1F_Map_Crop = T1F_Map; 
T1F_Map_Crop(~ROI) = 0;
figure(3); subplot(2,3,1); imagesc(flipud(T1F_Map_Crop)); cb = colorbar; cb.Label.String = '(s)'; set(cb,'FontSize',20); axis off; colormap(inferno); tt = title('T_{1F}'); tt.FontSize = 18; pbaspect([1.11 1.38 1])

T1S_Map = vec2mat(T1S_MapVector,length(x_vector));
T1S_Map_Crop = T1S_Map; 
T1S_Map_Crop(~ROI) = 0;
figure(3); subplot(2,3,2); imagesc(flipud(T1S_Map_Crop)); cb = colorbar; cb.Label.String = '(s)'; set(cb,'FontSize',20); axis off; colormap(inferno); tt = title('T_{1S}'); tt.FontSize = 18; pbaspect([1.11 1.38 1])

T2F_Map = vec2mat(T2F_MapVector,length(x_vector));
T2F_Map_Crop = T2F_Map; 
T2F_Map_Crop(~ROI) = 0;
figure(3); subplot(2,3,3); imagesc(flipud(T2F_Map_Crop)); cb = colorbar; cb.Label.String = '(s)'; set(cb,'FontSize',20); axis off; colormap(inferno); tt = title('T_{2F}'); tt.FontSize = 18; pbaspect([1.11 1.38 1])

T2S_Map = vec2mat(T2S_MapVector,length(x_vector));
T2S_Map_Crop = T2S_Map; 
T2S_Map_Crop(~ROI) = 0;
figure(3); subplot(2,3,4); imagesc(flipud(T2S_Map_Crop)); cb = colorbar; cb.Label.String = '(s)'; set(cb,'FontSize',20); axis off; colormap(inferno); tt = title('T_{2S}'); tt.FontSize = 18; pbaspect([1.11 1.38 1])

M0F_Map_Crop = M0F_Map; 
M0F_Map_Crop(~ROI) = 0;
figure(3); subplot(2,3,5); imagesc(flipud(M0F_Map_Crop)); cb = colorbar; set(cb,'FontSize',20); axis off; colormap(inferno); tt = title('MWF'); tt.FontSize = 20; pbaspect([1.11 1.38 1])

kFS_Map = vec2mat(kFS_MapVector,length(x_vector));
kFS_Map_Crop = kFS_Map; 
kFS_Map_Crop(~ROI) = 0;
figure(3); subplot(2,3,6); imagesc(flipud(kFS_Map_Crop)); cb = colorbar; cb.Label.String = '(s^{-1})'; set(cb,'FontSize',20); axis off; colormap(inferno); tt = title('k_{FS}'); tt.FontSize = 18; pbaspect([1.11 1.38 1])