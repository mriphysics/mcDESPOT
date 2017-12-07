function lnP = logpdf_2P (x, Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma)

if any(x < Lower) || any(x > Upper)
    lnP = -inf;
    return
end

if x(1) < x(2)
    lnP = -inf;
    return
end

if x(5) < x(6)
    lnP = -inf;
    return
end

Sig = [SPGR_steady_state_M0(FA_SPGR,TR_SPGR,'T1_S',x(1),'T1_F',x(2),'M0_F',x(3),'k_FS',x(4)); ...
    SSFP_steady_state_M0(FA_SSFP,TR_SSFP,'T1_S',x(1),'T2_S',x(5),'T1_F',x(2),'T2_F',x(6),'M0_F',x(3),'k_FS',x(4)); ...
    SSFP_steady_state_180_M0(FA_SSFP,TR_SSFP,'T1_S',x(1),'T2_S',x(5),'T1_F',x(2),'T2_F',x(6),'M0_F',x(3),'k_FS',x(4))];

lnP = -(sum(Sig - Data).^2)./(2*Sigma^2);

end