%%% Generate signal candidates for lsqnonlin optimisation. %%%

function Candidates = LM_Sigs(x,TR_SPGR,TR_SSFP, FA_SPGR, FA_SSFP0, FA_SSFP180)

[Mss_SPGR, Mss_SSFP0, Mss_SSFP180] = MEX_mcDESPOT_B0_DE_DiffFAs(x(1), x(2), x(3), x(4), x(5), x(6), TR_SPGR, TR_SSFP, x(7), FA_SPGR, FA_SSFP0, FA_SSFP180);

SSFP_Signals = [Mss_SSFP180, Mss_SSFP0];

Candidates = [Mss_SPGR./mean(Mss_SPGR,2) , SSFP_Signals./mean(SSFP_Signals,2)].';
    
end