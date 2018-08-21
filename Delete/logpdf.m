%% Cost-function calculation. Specify model and phase-cycling pattern.

function lnP = logpdf(x, Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma)

if any(x < Lower) || any(x > Upper)
    lnP = -inf;
    return
end

% Does not assume dynamic equilibrium, so fits for both M0s and both ks.
[Mss_SPGR, Mss_SSFP0, Mss_SSFP180] = MEX_mcDESPOT_B0_DiffSSFPFAs(x(2), x(1), x(8), x(7), x(5), x(6), x(3), x(4), TR_SPGR, TR_SSFP, x(9), rad2deg(FA_SPGR), rad2deg(FA_SSFP0), rad2deg(FA_SSFP180));
% Does assume dynamic equilibrium, so fits for kFS and M0F only.
%[Mss_SPGR, Mss_SSFP, Mss_SSFP180] = MEX_mcDESPOT_B0_DE_DiffFAs(x(2), x(1), x(8), x(7), x(5), x(3), TR_SPGR, TR_SSFP, x(9), rad2deg(FA_SPGR), rad2deg(FA_SSFP0), rad2deg(FA_SSFP180));
Mss_SSFP = [Mss_SSFP0 , Mss_SSFP180];

Data_SPGR = Data(1:numel(Mss_SPGR));
Data_SSFP = Data(numel(Mss_SPGR)+1:end);

Data_SPGR = (Data_SPGR * mean(Mss_SPGR))./mean(Data_SPGR);
Data_SSFP = (Data_SSFP * mean(Mss_SSFP))./mean(Data_SSFP);

Sig = [Mss_SPGR, Mss_SSFP]';
Data = [Data_SPGR; Data_SSFP];

lnP = -(sum((Sig - Data).^2))./(2*Sigma^2);

end

% PC1 = 0; PC2 = pi; PC3 = pi/2; PC4 = 3*(pi/2)? 
%[~, Mss_SSFP90, ~] = MEX_mcDESPOT_B0_DiffSSFPFAs(x(2), x(1), x(8), x(7), x(5), x(6), x(3), x(4), TR_SPGR, TR_SSFP, pi/2, rad2deg(FA_SPGR), rad2deg(FA_SSFP90), rad2deg(FA_SSFP180));
%[~, Mss_SSFP270, ~] = MEX_mcDESPOT_B0_DiffSSFPFAs(x(2), x(1), x(8), x(7), x(5), x(6), x(3), x(4), TR_SPGR, TR_SSFP, (3*pi/2), rad2deg(FA_SPGR), rad2deg(FA_SSFP270), rad2deg(FA_SSFP180));
