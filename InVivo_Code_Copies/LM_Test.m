%%% Fit mcDESPOT model to in vivo data using lsqnonlin. %%%

delete(gcp('nocreate')); c = parcluster('local'); c.NumWorkers = 8; parpool(c, c.NumWorkers);

Params_Est = zeros(length(Fitted_Coords),7);

lb = [0.2,0.8,0.002,0.06,0.5,0,0]; ub = [0.7,2.0,0.040,0.16,20,0.5,2*pi];

x0 = [((lb(1)+ub(1))/2),((lb(2)+ub(2))/2),((lb(3)+ub(3))/2),((lb(4)+ub(4))/2),((lb(5)+ub(5))/2),((lb(6)+ub(6))/2),((lb(7)+ub(7))/2)];

options = optimset('lsqnonlin'); options.Display = 'off'; options.TolFun = 1e-12; options.TolX = 1e-12; options.MaxIter = 10000; options.Algoritm = 'trust-region-reflective';

Fitted_FA_SPGR(2:36,:) = 1e-10; Fitted_FA_SSFP0(2:36,:) = 1e-10; Fitted_FA_SSFP180(2:36,:) = 1e-10;

parfor pix = 1:length(Fitted_Coords)
    
    Function = @(x)((LM_Sigs(x,TR_SPGR,TR_SSFP, Fitted_FA_SPGR(pix,:), Fitted_FA_SSFP0(pix,:), Fitted_FA_SSFP180(pix,:))) - Fitted_Data(:,pix));
    
    Params_Est(pix,:) = lsqnonlin(Function,x0,lb,ub,options);
    
end

%%% Run this to get data in correct format for Figure 10.
MWF_LM = zeros(size(Coords,1),1);
for pp = 1:length(Indices_Fitted)
   
    MWF_LM(Indices_Fitted(pp),1) = Params_Est(pp,6);

end

Map_BouhraraLSQ = flipud(vec2mat(MWF_LM,length(x_vector)));
%%%

ROI = roipoly(Map_BouhraraLSQ); 
MapCrop_BouhraraLSQ = Map_BouhraraLSQ; MapCrop_BouhraraLSQ(~ROI) = 0;
imagesc(MapCrop_BouhraraLSQ); cb1 = colorbar; axis off; colormap(magma); tt = title('lsqnonlin'); tt.FontSize = 16; pbaspect([1.11 1.38 1]); caxis([0 0.35]); cb1.FontSize = 12;
