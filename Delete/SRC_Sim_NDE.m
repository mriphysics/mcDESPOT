% SRC for pre-corrected FAs and pre-concatenated signals.

function [T1S_Sol, T1F_Sol, M0F_Sol, M0S_Sol, kFS_Sol, kSF_Sol, T2S_Sol, T2F_Sol] = SRC_Sim_NDE(Trials, Iterations, N, Runs, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data)

% Initialise loop variables.
T1S_SRC = zeros(Runs,1); T1F_SRC = zeros(Runs,1); M0F_SRC = zeros(Runs,1); M0S_SRC = zeros(Runs,1); T2S_SRC = zeros(Runs,1); T2F_SRC = zeros(Runs,1); kFS_SRC = zeros(Runs,1); kSF_SRC = zeros(Runs,1);
Candidates = zeros((length(FA_SPGR) + length(FA_SSFP0) + length(FA_SSFP180)),Trials);
Mss_SPGR = zeros(length(FA_SPGR),Trials); Mss_SSFP = zeros(length(FA_SSFP0),Trials); Mss_SSFP180 = zeros(length(FA_SSFP180),Trials);
SRSQ = zeros(Trials,1);
T1S_Sub = zeros(N,1); T1F_Sub = zeros(N,1); T2S_Sub = zeros(N,1); T2F_Sub = zeros(N,1); M0F_Sub = zeros(N,1); M0S_Sub = zeros(N,1); kFS_Sub = zeros(N,1); kSF_Sub = zeros(N,1);

T1S_LB = zeros(Iterations+1,1); T1S_UB = zeros(Iterations+1,1);
T2S_LB = zeros(Iterations+1,1); T2S_UB = zeros(Iterations+1,1);
T1F_LB = zeros(Iterations+1,1); T1F_UB = zeros(Iterations+1,1);
T2F_LB = zeros(Iterations+1,1); T2F_UB = zeros(Iterations+1,1);
M0F_LB = zeros(Iterations+1,1); M0F_UB = zeros(Iterations+1,1);
M0S_LB = zeros(Iterations+1,1); M0S_UB = zeros(Iterations+1,1);
kFS_LB = zeros(Iterations+1,1); kFS_UB = zeros(Iterations+1,1);
kSF_LB = zeros(Iterations+1,1); kSF_UB = zeros(Iterations+1,1);

for nn = 1:Runs
    
    T1S_LB(1) = 0.25; T1S_UB(1) = 1.5; T2S_LB(1) = 0.04; T2S_UB(1) = 0.15;
    T1F_LB(1) = 0.1; T1F_UB(1) = 0.8; T2F_LB(1) = 0.01; T2F_UB(1) = 0.03;
    M0F_LB(1) = 0.001; M0F_UB(1) = 0.35; M0S_LB(1) = 0.001; M0S_UB(1) = 0.8;
    kFS_LB(1) = (1/0.6); kFS_UB(1) = 40; kSF_LB(1) = (1/0.6); kSF_UB(1) = 20; % Wood.

    tic
    for it = 1:Iterations
        disp(['Iteration ', num2str(it)])
        
        if it == 1
            
            % Randomly choose parameter values from uniform distribution.
            T1S_Pick = (T1S_UB(it) - T1S_LB(it)) .* rand(Trials,1) + T1S_LB(it); T2S_Pick = (T2S_UB(it) - T2S_LB(it)) .* rand(Trials,1) + T2S_LB(it);
            T1F_Pick = (T1F_UB(it) - T1F_LB(it)) .* rand(Trials,1) + T1F_LB(it); T2F_Pick = (T2F_UB(it) - T2F_LB(it)) .* rand(Trials,1) + T2F_LB(it);
            M0F_Pick = (M0F_UB(it) - M0F_LB(it)) .* rand(Trials,1) + M0F_LB(it); kFS_Pick = (kFS_UB(it) - kFS_LB(it)) .* rand(Trials,1) + kFS_LB(it);
            M0S_Pick = (M0S_UB(it) - M0S_LB(it)) .* rand(Trials,1) + M0S_LB(it); kSF_Pick = (kSF_UB(it) - kSF_LB(it)) .* rand(Trials,1) + kSF_LB(it);
            
        else
            
            T1F_Pick = mean(T1F_Sub) + (std(T1F_Sub) .* trandn((((T1F_LB(it)*ones(Trials,1))-mean(T1F_Sub))./std(T1F_Sub)), (((T1F_UB(it)*ones(Trials,1))-mean(T1F_Sub))./std(T1F_Sub))));
            T1S_Pick = mean(T1S_Sub) + (std(T1S_Sub) .* trandn((((T1S_LB(it)*ones(Trials,1))-mean(T1S_Sub))./std(T1S_Sub)), (((T1S_UB(it)*ones(Trials,1))-mean(T1S_Sub))./std(T1S_Sub))));
            T2F_Pick = mean(T2F_Sub) + (std(T2F_Sub) .* trandn((((T2F_LB(it)*ones(Trials,1))-mean(T2F_Sub))./std(T2F_Sub)), (((T2F_UB(it)*ones(Trials,1))-mean(T2F_Sub))./std(T2F_Sub))));
            T2S_Pick = mean(T2S_Sub) + (std(T2S_Sub) .* trandn((((T2S_LB(it)*ones(Trials,1))-mean(T2S_Sub))./std(T2S_Sub)), (((T2S_UB(it)*ones(Trials,1))-mean(T2S_Sub))./std(T2S_Sub))));
            M0F_Pick = mean(M0F_Sub) + (std(M0F_Sub) .* trandn((((M0F_LB(it)*ones(Trials,1))-mean(M0F_Sub))./std(M0F_Sub)), (((M0F_UB(it)*ones(Trials,1))-mean(M0F_Sub))./std(M0F_Sub))));
            kFS_Pick = mean(kFS_Sub) + (std(kFS_Sub) .* trandn((((kFS_LB(it)*ones(Trials,1))-mean(kFS_Sub))./std(kFS_Sub)), (((kFS_UB(it)*ones(Trials,1))-mean(kFS_Sub))./std(kFS_Sub))));
            M0S_Pick = mean(M0S_Sub) + (std(M0S_Sub) .* trandn((((M0S_LB(it)*ones(Trials,1))-mean(M0S_Sub))./std(M0S_Sub)), (((M0S_UB(it)*ones(Trials,1))-mean(M0S_Sub))./std(M0S_Sub))));
            kSF_Pick = mean(kSF_Sub) + (std(kSF_Sub) .* trandn((((kSF_LB(it)*ones(Trials,1))-mean(kSF_Sub))./std(kSF_Sub)), (((kSF_UB(it)*ones(Trials,1))-mean(kSF_Sub))./std(kSF_Sub))));    
            
        end
        
        for qq = 1:Trials
            
            % Generate a number of theoretical signals, normalise and concatenate.
            Mss_SPGR(:,qq) = SPGR_SteadyState(FA_SPGR,TR_SPGR,'T1_S',T1S_Pick(qq),'T1_F',T1F_Pick(qq),'M0_F',M0F_Pick(qq),'M0_S',M0S_Pick(qq),'k_FS',kFS_Pick(qq),'k_SF',kSF_Pick(qq));
            Mss_SSFP(:,qq) = SSFP_SteadyState(FA_SSFP0,TR_SSFP,0,'T1_S',T1S_Pick(qq),'T2_S',T2S_Pick(qq),'T1_F',T1F_Pick(qq),'T2_F',T2F_Pick(qq),'M0_F',M0F_Pick(qq),'M0_S',M0S_Pick(qq),'k_FS',kFS_Pick(qq),'k_SF',kSF_Pick(qq));
            Mss_SSFP180(:,qq) = SSFP_SteadyState(FA_SSFP180,TR_SSFP,pi,'T1_S',T1S_Pick(qq),'T2_S',T2S_Pick(qq),'T1_F',T1F_Pick(qq),'T2_F',T2F_Pick(qq),'M0_F',M0F_Pick(qq),'M0_S',M0S_Pick(qq),'k_FS',kFS_Pick(qq),'k_SF',kSF_Pick(qq));
            Candidates(:,qq) = [Mss_SPGR(:,qq) ; Mss_SSFP(:,qq) ; Mss_SSFP180(:,qq)];
            SRSQ(qq) = norm(Candidates(:,qq) - Data).^2;

            % Does not assume dynamic equilibrium, so fits for both M0s and both ks.
            %[Mss_SPGR, Mss_SSFP, Mss_SSFP180] = MEX_mcDESPOT_B0_DiffSSFPFAs(x(2), x(1), x(8), x(7), x(5), x(6), x(3), x(4), TR_SPGR, TR_SSFP, x(9), rad2deg(FA_SPGR), rad2deg(FA_SSFP0), rad2deg(FA_SSFP180));
            % Does assume dynamic equilibrium, so fits for kFS and M0F only.
            %[Mss_SPGR, Mss_SSFP, Mss_SSFP180] = MEX_mcDESPOT_B0_DE_DiffFAs(x(2), x(1), x(8), x(7), x(5), x(3), TR_SPGR, TR_SSFP, x(9), rad2deg(FA_SPGR), rad2deg(FA_SSFP0), rad2deg(FA_SSFP180));
            
        end
        
        % Sort residuals in ascending order to determine indices of N-lowest residual combinations.
        [~, Index] = sort(SRSQ);
        
        % Locate indices in existing parameter sets.
        for pp = 1:N
            
            T1S_Sub(pp) = T1S_Pick(Index(pp)); T1F_Sub(pp) = T1F_Pick(Index(pp));
            T2S_Sub(pp) = T2S_Pick(Index(pp)); T2F_Sub(pp) = T2F_Pick(Index(pp));
            M0F_Sub(pp) = M0F_Pick(Index(pp)); kFS_Sub(pp) = kFS_Pick(Index(pp));
            M0S_Sub(pp) = M0S_Pick(Index(pp)); kSF_Sub(pp) = kSF_Pick(Index(pp));
            
        end
        
        % Extract new search-space bounds from each subset.
        T1S_LB(it+1) = min(T1S_Sub) ; T1S_UB(it+1) = max(T1S_Sub);
        T1F_LB(it+1) = min(T1F_Sub) ; T1F_UB(it+1) = max(T1F_Sub);
        M0F_LB(it+1) = min(M0F_Sub) ; M0F_UB(it+1) = max(M0F_Sub);
        M0S_LB(it+1) = min(M0S_Sub) ; M0S_UB(it+1) = max(M0S_Sub);
        T2S_LB(it+1) = min(T2S_Sub) ; T2S_UB(it+1) = max(T2S_Sub);
        T2F_LB(it+1) = min(T2F_Sub) ; T2F_UB(it+1) = max(T2F_Sub);
        kFS_LB(it+1) = min(kFS_Sub) ; kFS_UB(it+1) = max(kFS_Sub);
        kSF_LB(it+1) = min(kSF_Sub) ; kSF_UB(it+1) = max(kSF_Sub);
        
        % Expand search-space slightly to avoid 'inadvertent over-contraction'.
        T1S_expn = (T1S_UB(it+1) - T1S_LB(it+1))/N; T1F_expn = (T1F_UB(it+1) - T1F_LB(it+1))/N;
        T2S_expn = (T2S_UB(it+1) - T2S_LB(it+1))/N; T2F_expn = (T2F_UB(it+1) - T2F_LB(it+1))/N;
        M0F_expn = (M0F_UB(it+1) - M0F_LB(it+1))/N; kFS_expn = (kFS_UB(it+1) - kFS_LB(it+1))/N;
        M0S_expn = (M0S_UB(it+1) - M0S_LB(it+1))/N; kSF_expn = (kSF_UB(it+1) - kSF_LB(it+1))/N;
        
        T1S_LB(it+1) = T1S_LB(it+1) - T1S_expn; T1S_UB(it+1) = T1S_UB(it+1) + T1S_expn;
        T1F_LB(it+1) = T1F_LB(it+1) - T1F_expn; T1F_UB(it+1) = T1F_UB(it+1) + T1F_expn;
        T2S_LB(it+1) = T2S_LB(it+1) - T2S_expn; T2S_UB(it+1) = T2S_UB(it+1) + T2S_expn;
        T2F_LB(it+1) = T2F_LB(it+1) - T2F_expn; T2F_UB(it+1) = T2F_UB(it+1) + T2F_expn;
        M0F_LB(it+1) = M0F_LB(it+1) - M0F_expn; M0F_UB(it+1) = M0F_UB(it+1) + M0F_expn;
        kFS_LB(it+1) = kFS_LB(it+1) - kFS_expn; kFS_UB(it+1) = kFS_UB(it+1) + kFS_expn;
        M0S_LB(it+1) = M0S_LB(it+1) - M0S_expn; M0S_UB(it+1) = M0S_UB(it+1) + M0S_expn;
        kSF_LB(it+1) = kSF_LB(it+1) - kSF_expn; kSF_UB(it+1) = kSF_UB(it+1) + kSF_expn;
        
        % Set lower bound to zero if less than zero after expansion.
        if T1S_LB(it+1) < 0
            T1S_LB(it+1) = 0;
        end
        if T2S_LB(it+1) < 0
            T2S_LB(it+1) = 0;
        end
        if T1F_LB(it+1) < 0
            T1F_LB(it+1) = 0;
        end
        if T2F_LB(it+1) < 0
            T2F_LB(it+1) = 0;
        end
        if M0F_LB(it+1) < 0
            M0F_LB(it+1) = 0;
        end
        if kFS_LB(it+1) < 0
            kFS_LB(it+1) = 0;
        end
        if M0S_LB(it+1) < 0
            M0S_LB(it+1) = 0;
        end
        if kSF_LB(it+1) < 0
            kSF_LB(it+1) = 0;
        end
        
        if (((max(T1S_Sub)-min(T1S_Sub))/max(T1S_Sub)) < 0.01 && ((max(T1F_Sub)-min(T1F_Sub))/max(T1F_Sub)) < 0.01 && ((max(M0F_Sub)-min(M0F_Sub))/max(M0F_Sub)) < 0.01 && ((max(kFS_Sub)-min(kFS_Sub))/max(kFS_Sub)) < 0.01 && ((max(T2S_Sub)-min(T2S_Sub))/max(T2S_Sub)) < 0.01 && ((max(T2F_Sub)-min(T2F_Sub))/max(T2F_Sub)) < 0.01 && ((max(M0S_Sub)-min(M0S_Sub))/max(M0S_Sub)) < 0.01 && ((max(kSF_Sub)-min(kSF_Sub))/max(kSF_Sub)) < 0.01)
            break
        end
        
    end
    toc
    % Return result as mean of top-five solutions.
    T1S_SRC(nn) = mean(T1S_Sub(1:5)); T1F_SRC(nn) = mean(T1F_Sub(1:5)); T2S_SRC(nn) = mean(T2S_Sub(1:5)); T2F_SRC(nn) = mean(T2F_Sub(1:5)); M0F_SRC(nn) = mean(M0F_Sub(1:5)); kFS_SRC(nn) = mean(kFS_Sub(1:5)); M0S_SRC(nn) = mean(M0S_Sub(1:5)); kSF_SRC(nn) = mean(kSF_Sub(1:5));
    
    % Second-stage optimisation.
    %Options = optimoptions(@fmincon, 'MaxFunctionEvaluations', 10000, 'MaxIterations', 10000, 'OptimalityTolerance', 5e-8, 'Display', 'none');
    %x0(nn,:) = [T1S_SRC(nn) T1F_SRC(nn) M0F_SRC(nn) kFS_SRC(nn) T2S_SRC(nn) T2F_SRC(nn)];
    %Sig_SPGR = @(x)(SPGR_SteadyState(FA_SPGR,TR_SPGR,'T1_S',x(1),'T1_F',x(2),'M0_F',x(3),'k_FS',x(4)));
    %Sig_SSFP0 = @(x)(SSFP_SteadyState(FA_SSFP0,TR_SSFP,0,'T1_S',x(1),'T2_S',x(5),'T1_F',x(2),'T2_F',x(6),'M0_F',x(3),'k_FS',x(4)));
    %Sig_SSFP180 = @(x)(SSFP_SteadyState(FA_SSFP180,TR_SSFP,pi,'T1_S',x(1),'T2_S',x(5),'T1_F',x(2),'T2_F',x(6),'M0_F',x(3),'k_FS',x(4)));
    %Sig = @(x)([Sig_SPGR(x) ; Sig_SSFP0(x) ; Sig_SSFP180(x)]);
    %CF = @(x) double(norm(Sig(x) - Data).^2);
    %[Sol(nn,:),Fval(nn),Exitflag(nn),~] = fmincon(CF,x0(nn,:),[],[],[],[],[T1S_LB(it+1) T1F_LB(it+1) M0F_LB(it+1) kFS_LB(it+1) T2S_LB(it+1) T2F_LB(it+1)],[T1S_UB(it+1) T1F_UB(it+1) M0F_UB(it+1) kFS_UB(it+1) T2S_UB(it+1) T2F_UB(it+1)],[],Options);
    %if Exitflag(nn) ~= 0
    %    T1S_app(nn) = Sol(nn,1); T1F_app(nn) = Sol(nn,2); M0F_app(nn) = Sol(nn,3); kFS_app(nn) = Sol(nn,4); T2S_app(nn) = Sol(nn,5); T2F_app(nn) = Sol(nn,6);
    %end
        
end

    T1S_Sol= mean(T1S_SRC); T1F_Sol = mean(T1F_SRC); T2S_Sol = mean(T2S_SRC); T2F_Sol = mean(T2F_SRC); M0F_Sol = mean(M0F_SRC); kFS_Sol = mean(kFS_SRC); M0S_Sol = mean(M0S_SRC); kSF_Sol = mean(kSF_SRC);

end