function lnP = logpdf (x, Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma)

if any(x < Lower) || any(x > Upper)
    lnP = -inf;
    return
end

if x(1) < x(2)
    lnP = -inf;
    return
end

if x(7) < x(8)
    lnP = -inf;
    return
end

Sig = [SPGR_steady_state_M0(FA_SPGR,TR_SPGR,'T1_S',x(1),'T1_F',x(2),'M0_F',x(3),'M0_S',x(4),'k_FS',x(5),'k_SF',x(6)); ...
    SSFP_steady_state_M0(FA_SSFP,TR_SSFP,'T1_S',x(1),'T2_S',x(7),'T1_F',x(2),'T2_F',x(8),'M0_F',x(3),'M0_S',x(4),'k_FS',x(5),'k_SF',x(6)); ...
    SSFP_steady_state_180_M0(FA_SSFP,TR_SSFP,'T1_S',x(1),'T2_S',x(7),'T1_F',x(2),'T2_F',x(8),'M0_F',x(3),'M0_S',x(4),'k_FS',x(5),'k_SF',x(6))];

lnP = -(sum(Sig - Data).^2)./(2*Sigma^2);

end