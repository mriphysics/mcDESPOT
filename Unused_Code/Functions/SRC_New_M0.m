%% mcDESPOT processing as per Deoni (2015) using steady-state equations from Liu (2016).

close all
clear all

Trials = 5000; Iterations = 7; N = 50; Runs = 1;

% Define sequence and ground-truth tissue parameters. Increased FA-range for high signal curve definition.
TR_SPGR = 5.2e-3; TR_SSFP = 5.2e-3;
FA_SPGR = deg2rad([3 4 5 6 7 9 13 18]); FA_SSFP = deg2rad([2 5 10 15 20 30 40 50]);
T1_S = 1.15; T2_S = 0.08; T1_F = 0.4; T2_F = 0.02; M0_F = 0.25; M0_S = 0.75; k_FS = 9; T1_B = 1; T2_B = 8e-6; M0_B = 0.2; k_SB = 5; k_FB = 5;

% Initialise loop variables.
T1S_final = zeros(length(M0_B),1); T1F_final = zeros(length(M0_B),1); T2S_final = zeros(length(M0_B),1); T2F_final = zeros(length(M0_B),1); M0F_final = zeros(length(M0_B),1); M0S_final = zeros(length(M0_B),1); kFS_final = zeros(length(M0_B),1); kSF_final = zeros(length(M0_B),1);
T1S_SRC = zeros(Runs,length(M0_B)); T1F_SRC = zeros(Runs,length(M0_B)); M0F_SRC = zeros(Runs,length(M0_B)); M0S_SRC = zeros(Runs,length(M0_B)); kFS_SRC = zeros(Runs,length(M0_B)); kSF_SRC = zeros(Runs,length(M0_B)); T2S_SRC = zeros(Runs,length(M0_B)); T2F_SRC = zeros(Runs,length(M0_B));
T1S_app = zeros(Runs,length(M0_B)); T1F_app = zeros(Runs,length(M0_B)); M0F_app = zeros(Runs,length(M0_B)); M0S_app = zeros(Runs,length(M0_B)); kFS_app = zeros(Runs,length(M0_B)); kSF_app = zeros(Runs,length(M0_B)); T2S_app = zeros(Runs,length(M0_B)); T2F_app = zeros(Runs,length(M0_B));

ErrorTP_T1S = zeros(length(M0_B),1); ErrorTP_T2S = zeros(length(M0_B),1);
ErrorTP_T1F = zeros(length(M0_B),1); ErrorTP_T2F = zeros(length(M0_B),1);
ErrorTP_M0F = zeros(length(M0_B),1); ErrorTP_kFS = zeros(length(M0_B),1);
PrecisionTP_T1S = zeros(length(M0_B),1); PrecisionTP_T2S = zeros(length(M0_B),1);
PrecisionTP_T1F = zeros(length(M0_B),1); PrecisionTP_T2F = zeros(length(M0_B),1);
PrecisionTP_M0F = zeros(length(M0_B),1); PrecisionTP_kFS = zeros(length(M0_B),1);

kSF_Prediction = zeros(length(M0_B),1); kFS_Prediction = zeros(length(M0_B),1);
R1F_Prediction = zeros(length(M0_B),1); R1S_Prediction = zeros(length(M0_B),1);
M0F_Prediction = zeros(length(M0_B),1); M0S_Prediction = zeros(length(M0_B),1);

kSF_Prediction2 = zeros(length(M0_B),1); kFS_Prediction2 = zeros(length(M0_B),1);
R1F_Prediction2 = zeros(length(M0_B),1); R1S_Prediction2 = zeros(length(M0_B),1);
M0F_Prediction2 = zeros(length(M0_B),1); M0S_Prediction2 = zeros(length(M0_B),1);

for ss = 1:length(M0_B)

    % Initiate loop variables.
    Candidates = zeros((length(FA_SPGR) + length(FA_SSFP) + length(FA_SSFP)),Trials);
    Mss_SPGR = zeros(length(FA_SPGR),Trials); Mss_SSFP = zeros(length(FA_SSFP),Trials); Mss_SSFP_180 = zeros(length(FA_SSFP),Trials);
    SRSQ = zeros(Trials,1);
    T1S_Sub = zeros(N,1); T1F_Sub = zeros(N,1); T2S_Sub = zeros(N,1); T2F_Sub = zeros(N,1); M0F_Sub = zeros(N,1); M0S_Sub = zeros(N,1); kFS_Sub = zeros(N,1); kSF_Sub = zeros(N,1);
    exitflag = zeros(Runs,1); fval = zeros(Runs,1); Sol = zeros(Runs,8); x0 = zeros(Runs,8);
    
    % Define ground-truth signals and normalise.
    SPGR_Data = CS_SPGR_steady_state_MT_M0(FA_SPGR,TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'T1_B',T1_B,'M0_B',M0_B(ss),'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
    SSFP_Data = CS_SSFP_steady_state_MT_M0(FA_SSFP,TR_SSFP,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B(ss),'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
    SSFP_Data_180 = CS_SSFP_steady_state_180_MT_M0(FA_SSFP,TR_SSFP,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B(ss),'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
    Data = [SPGR_Data ; SSFP_Data ; SSFP_Data_180];
    
    % Calculation of noise factor. Specify desired SNR.
    % M0 = 1; SNR = 100; NF = (0.1 * M0)/SNR;
    
    T1S_LB = zeros(Iterations+1,1); T1S_UB = zeros(Iterations+1,1);
    T2S_LB = zeros(Iterations+1,1); T2S_UB = zeros(Iterations+1,1);
    T1F_LB = zeros(Iterations+1,1); T1F_UB = zeros(Iterations+1,1);
    T2F_LB = zeros(Iterations+1,1); T2F_UB = zeros(Iterations+1,1);
    M0F_LB = zeros(Iterations+1,1); M0F_UB = zeros(Iterations+1,1);
    M0S_LB = zeros(Iterations+1,1); M0S_UB = zeros(Iterations+1,1);
    kFS_LB = zeros(Iterations+1,1); kFS_UB = zeros(Iterations+1,1);
    kSF_LB = zeros(Iterations+1,1); kSF_UB = zeros(Iterations+1,1);
    
    for nn = 1:Runs
        
        disp(['Run Number: ', num2str(nn)])
        
        % Specify initial upper and lower bounds.
        T1S_LB(1) = 0.3; T1S_UB(1) = 1.5; T2S_LB(1) = 0.04; T2S_UB(1) = 0.15;
        T1F_LB(1) = 0.1; T1F_UB(1) = 0.8; T2F_LB(1) = 0.01; T2F_UB(1) = 0.03;
        M0F_LB(1) = 0.001; M0F_UB(1) = 0.35; M0S_LB(1) = 0.3; M0S_UB(1) = 0.8; 
        kFS_LB(1) = (1/0.6); kFS_UB(1) = 40; kSF_LB(1) = (1/0.6); kSF_UB(1) = 40;
        
        % SPGR_Data_Noisy = SPGR_Data + complex(NF * randn(size(SPGR_Data)), NF * randn(size(SPGR_Data)));
        % SSFP_Data_Noisy = SSFP_Data + complex(NF * randn(size(SSFP_Data)), NF * randn(size(SSFP_Data)));
        % SSFP_Data_Noisy_180 = SSFP_Data_180 + complex(NF * randn(size(SSFP_Data_180)), NF * randn(size(SSFP_Data_180)));
        % Concatenate signals for mcDESPOT processing. Choose noisy or noiseless.
        % Data = [SPGR_Data_Noisy ; SSFP_Data_Noisy ; SSFP_Data_Noisy_180];
        
        tic
        for jj = 1:Iterations
            
            if jj == 1
                
                % Randomly choose parameter values from uniform distribution.
                T1S_Pick = (T1S_UB(jj) - T1S_LB(jj)) .* rand(Trials,1) + T1S_LB(jj); T2S_Pick = (T2S_UB(jj) - T2S_LB(jj)) .* rand(Trials,1) + T2S_LB(jj);
                T1F_Pick = (T1F_UB(jj) - T1F_LB(jj)) .* rand(Trials,1) + T1F_LB(jj); T2F_Pick = (T2F_UB(jj) - T2F_LB(jj)) .* rand(Trials,1) + T2F_LB(jj);
                M0F_Pick = (M0F_UB(jj) - M0F_LB(jj)) .* rand(Trials,1) + M0F_LB(jj); M0S_Pick = (M0S_UB(jj) - M0S_LB(jj)) .* rand(Trials,1) + M0S_LB(jj);
                kFS_Pick = (kFS_UB(jj) - kFS_LB(jj)) .* rand(Trials,1) + kFS_LB(jj); kSF_Pick = (kSF_UB(jj) - kSF_LB(jj)) .* rand(Trials,1) + kSF_LB(jj);
                
            else
                
                NormDist_T1S = makedist('normal', 'mu', mean(T1S_Sub),'sigma',std(T1S_Sub));
                NormDist_T1S = truncate(NormDist_T1S, T1S_LB(jj), T1S_UB(jj));
                NormDist_T2S = makedist('normal', 'mu', mean(T2S_Sub),'sigma',std(T2S_Sub));
                NormDist_T2S = truncate(NormDist_T2S, T2S_LB(jj), T2S_UB(jj));
                NormDist_T1F = makedist('normal', 'mu', mean(T1F_Sub),'sigma',std(T1F_Sub));
                NormDist_T1F = truncate(NormDist_T1F, T1F_LB(jj), T1F_UB(jj));
                NormDist_T2F = makedist('normal', 'mu', mean(T2F_Sub),'sigma',std(T2F_Sub));
                NormDist_T2F = truncate(NormDist_T2F, T2F_LB(jj), T2F_UB(jj));
                NormDist_M0F = makedist('normal', 'mu', mean(M0F_Sub),'sigma',std(M0F_Sub));
                NormDist_M0F = truncate(NormDist_M0F, M0F_LB(jj), M0F_UB(jj));
                NormDist_M0S = makedist('normal', 'mu', mean(M0S_Sub),'sigma',std(M0S_Sub));
                NormDist_M0S = truncate(NormDist_M0S, M0S_LB(jj), M0S_UB(jj));
                NormDist_kFS = makedist('normal', 'mu', mean(kFS_Sub),'sigma',std(kFS_Sub));
                NormDist_kFS = truncate(NormDist_kFS, kFS_LB(jj), kFS_UB(jj));
                NormDist_kSF = makedist('normal', 'mu', mean(kSF_Sub),'sigma',std(kSF_Sub));
                NormDist_kSF = truncate(NormDist_kSF, kSF_LB(jj), kSF_UB(jj));
                
                for tt = 1:Trials
                    T1S_Pick(tt) = random(NormDist_T1S);
                    T2S_Pick(tt) = random(NormDist_T2S);
                    T1F_Pick(tt) = random(NormDist_T1F);
                    T2F_Pick(tt) = random(NormDist_T2F);
                    M0F_Pick(tt) = random(NormDist_M0F);
                    M0S_Pick(tt) = random(NormDist_M0S);
                    kFS_Pick(tt) = random(NormDist_kFS);
                    kSF_Pick(tt) = random(NormDist_kSF);
                end
                
                
            end
            
            for ii = 1:Trials

                % Generate a number of theoretical signals, normalise and concatenate.
                % Non-MT.   
                Mss_SPGR(:,ii) = SPGR_steady_state_M0(FA_SPGR,TR_SPGR,'T1_S',T1S_Pick(ii),'T1_F',T1F_Pick(ii),'M0_F',M0F_Pick(ii),'M0_S',M0S_Pick(ii),'k_FS',kFS_Pick(ii),'k_SF',kSF_Pick(ii));
                Mss_SSFP(:,ii) = SSFP_steady_state_M0(FA_SSFP,TR_SSFP,'T1_S',T1S_Pick(ii),'T2_S',T2S_Pick(ii),'T1_F',T1F_Pick(ii),'T2_F',T2F_Pick(ii),'M0_F',M0F_Pick(ii),'M0_S',M0S_Pick(ii),'k_FS',kFS_Pick(ii),'k_SF',kSF_Pick(ii));
                Mss_SSFP_180(:,ii) = SSFP_steady_state_180_M0(FA_SSFP,TR_SSFP,'T1_S',T1S_Pick(ii),'T2_S',T2S_Pick(ii),'T1_F',T1F_Pick(ii),'T2_F',T2F_Pick(ii),'M0_F',M0F_Pick(ii),'M0_S',M0S_Pick(ii),'k_FS',kFS_Pick(ii),'k_SF',kSF_Pick(ii));
                Candidates(:,ii) = [Mss_SPGR(:,ii) ; Mss_SSFP(:,ii) ; Mss_SSFP_180(:,ii)];
                % Calculate sum-of-squares residuals between each combination and acquired data.
                SRSQ(ii) = norm(Candidates(:,ii) - Data)^2;
                
            end
            
            % Sort residuals in ascending order to determine indices of N-lowest residual combinations.
            [Value, Index] = sort(SRSQ);
            
            % if jj <  Iterations
            
            % Locate indices in existing parameter sets.
            for ll = 1:N
                
                T1S_Sub(ll) = T1S_Pick(Index(ll)); T1F_Sub(ll) = T1F_Pick(Index(ll));
                T2S_Sub(ll) = T2S_Pick(Index(ll)); T2F_Sub(ll) = T2F_Pick(Index(ll));
                M0F_Sub(ll) = M0F_Pick(Index(ll)); M0S_Sub(ll) = M0S_Pick(Index(ll));
                kFS_Sub(ll) = kFS_Pick(Index(ll)); kSF_Sub(ll) = kSF_Pick(Index(ll));
                
            end
            
            % Extract new search-space bounds from each subset.
            T1S_LB(jj+1) = min(T1S_Sub) ; T1S_UB(jj+1) = max(T1S_Sub);
            T1F_LB(jj+1) = min(T1F_Sub) ; T1F_UB(jj+1) = max(T1F_Sub);
            M0F_LB(jj+1) = min(M0F_Sub) ; M0F_UB(jj+1) = max(M0F_Sub);
            M0S_LB(jj+1) = min(M0S_Sub) ; M0S_UB(jj+1) = max(M0S_Sub);
            kFS_LB(jj+1) = min(kFS_Sub) ; kFS_UB(jj+1) = max(kFS_Sub);
            kSF_LB(jj+1) = min(kSF_Sub) ; kSF_UB(jj+1) = max(kSF_Sub);
            T2S_LB(jj+1) = min(T2S_Sub) ; T2S_UB(jj+1) = max(T2S_Sub);
            T2F_LB(jj+1) = min(T2F_Sub) ; T2F_UB(jj+1) = max(T2F_Sub);
            
            % Expand search-space slightly to avoid 'inadvertent over-contraction'.
            T1S_expn = (T1S_UB(jj+1) - T1S_LB(jj+1))/N; T1F_expn = (T1F_UB(jj+1) - T1F_LB(jj+1))/N;
            T2S_expn = (T2S_UB(jj+1) - T2S_LB(jj+1))/N; T2F_expn = (T2F_UB(jj+1) - T2F_LB(jj+1))/N;
            M0F_expn = (M0F_UB(jj+1) - M0F_LB(jj+1))/N; M0S_expn = (M0S_UB(jj+1) - M0S_LB(jj+1))/N; 
            kFS_expn = (kFS_UB(jj+1) - kFS_LB(jj+1))/N; kSF_expn = (kSF_UB(jj+1) - kSF_LB(jj+1))/N;
            
            T1S_LB(jj+1) = T1S_LB(jj+1) - T1S_expn; T1S_UB(jj+1) = T1S_UB(jj+1) + T1S_expn;
            T1F_LB(jj+1) = T1F_LB(jj+1) - T1F_expn; T1F_UB(jj+1) = T1F_UB(jj+1) + T1F_expn;
            T2S_LB(jj+1) = T2S_LB(jj+1) - T2S_expn; T2S_UB(jj+1) = T2S_UB(jj+1) + T2S_expn;
            T2F_LB(jj+1) = T2F_LB(jj+1) - T2F_expn; T2F_UB(jj+1) = T2F_UB(jj+1) + T2F_expn;
            M0F_LB(jj+1) = M0F_LB(jj+1) - M0F_expn; M0F_UB(jj+1) = M0F_UB(jj+1) + M0F_expn;
            M0S_LB(jj+1) = M0S_LB(jj+1) - M0S_expn; M0S_UB(jj+1) = M0S_UB(jj+1) + M0S_expn;
            kFS_LB(jj+1) = kFS_LB(jj+1) - kFS_expn; kFS_UB(jj+1) = kFS_UB(jj+1) + kFS_expn;
            kSF_LB(jj+1) = kSF_LB(jj+1) - kSF_expn; kSF_UB(jj+1) = kSF_UB(jj+1) + kSF_expn;
            
            % Set lower bound to zero if less than zero after expansion.
            if T1S_LB(jj+1) < 0
                T1S_LB(jj+1) = 0;
            end
            if T2S_LB(jj+1) < 0
                T2S_LB(jj+1) = 0;
            end
            if T1F_LB(jj+1) < 0
                T1F_LB(jj+1) = 0;
            end
            if T2F_LB(jj+1) < 0
                T2F_LB(jj+1) = 0;
            end
            if M0F_LB(jj+1) < 0
                M0F_LB(jj+1) = 0;
            end
            if M0S_LB(jj+1) < 0
                M0S_LB(jj+1) = 0;
            end
            if kFS_LB(jj+1) < 0
                kFS_LB(jj+1) = 0;
            end
            if kSF_LB(jj+1) < 0
                kSF_LB(jj+1) = 0;
            end
            
            % end
            
            % Exit loop if difference between min and max of each parameter < 1%.
            if (((max(T1S_Sub)-min(T1S_Sub))/max(T1S_Sub)) < 0.01 && ((max(T1F_Sub)-min(T1F_Sub))/max(T1F_Sub)) < 0.01 && ((max(M0F_Sub)-min(M0F_Sub))/max(M0F_Sub)) < 0.01 && ((max(M0S_Sub)-min(M0S_Sub))/max(M0S_Sub)) < 0.01 && ((max(kFS_Sub)-min(kFS_Sub))/max(kFS_Sub)) < 0.01 && ((max(kSF_Sub)-min(kSF_Sub))/max(kSF_Sub)) < 0.01 && ((max(T2S_Sub)-min(T2S_Sub))/max(T2S_Sub)) < 0.01 && ((max(T2F_Sub)-min(T2F_Sub))/max(T2F_Sub)) < 0.01)
                break
            end
            
        end
        toc
        
        % Return result as mean of top-five solutions.
        T1S_SRC(nn,ss) = mean(T1S_Sub(1:5)); T1F_SRC(nn,ss) = mean(T1F_Sub(1:5)); T2S_SRC(nn,ss) = mean(T2S_Sub(1:5)); T2F_SRC(nn,ss) = mean(T2F_Sub(1:5)); M0F_SRC(nn,ss) = mean(M0F_Sub(1:5)); M0S_SRC(nn,ss) = mean(M0S_Sub(1:5)); kFS_SRC(nn,ss) = mean(kFS_Sub(1:5)); kSF_SRC(nn,ss) = mean(kSF_Sub(1:5));
        
        % OPTIMIZATION STEP II: Use guess as initial guess for local optim (fminsearch).
        tic
        % options = optimset('MaxFunEvals',1e10);
        x0(nn,:) = [T1S_SRC(nn,ss) T1F_SRC(nn,ss) M0F_SRC(nn,ss) M0S_SRC(nn,ss) kFS_SRC(nn,ss) kSF_SRC(nn,ss) T2S_SRC(nn,ss) T2F_SRC(nn,ss)];
        % Non-MT.
        Sig_SPGR = @(x)(SPGR_steady_state_M0(FA_SPGR,TR_SPGR,'T1_S',x(1),'T1_F',x(2),'M0_F',x(3),'M0_S',x(4),'k_FS',x(5),'k_SF',x(6)));
        Sig_SSFP = @(x)(SSFP_steady_state_M0(FA_SSFP,TR_SSFP,'T1_S',x(1),'T2_S',x(7),'T1_F',x(2),'T2_F',x(8),'M0_F',x(3),'M0_S',x(4),'k_FS',x(5),'k_SF',x(6)));
        Sig_SSFP180 = @(x)(SSFP_steady_state_180_M0(FA_SSFP,TR_SSFP,'T1_S',x(1),'T2_S',x(7),'T1_F',x(2),'T2_F',x(8),'M0_F',x(3),'M0_S',x(4),'k_FS',x(5),'k_SF',x(6)));
        Sig = @(x)([Sig_SPGR(x) ; Sig_SSFP(x) ; Sig_SSFP180(x)]);
        % Objective function should be same as that for SRC.
        CF = @(x) norm(Sig(x) - Data)^2;
        %[SolTP(nn,:),fval(nn),exitflag(nn),output] = fminsearch(cfTP,x0TP(nn,:));
        [Sol(nn,:),fval(nn),exitflag(nn),output] = fmincon(CF,x0(nn,:),[],[],[],[],[T1S_LB(jj+1) T1F_LB(jj+1) M0F_LB(jj+1) M0S_LB(jj+1) kFS_LB(jj+1) kSF_LB(jj+1) T2S_LB(jj+1) T2F_LB(jj+1)],[T1S_UB(jj+1) T1F_UB(jj+1) M0F_UB(jj+1) M0S_UB(jj+1) kFS_UB(jj+1) kSF_UB(jj+1) T2S_UB(jj+1) T2F_UB(jj+1)]);
        % Use an if-loop to remove non-convergent solutions?
        %if exitflag(nn) == 1;
        T1S_app(ii,ss) = Sol(1); T1F_app(ii,ss) = Sol(2); M0F_app(ii,ss) = Sol(3); M0S_app(ii,ss) = Sol(4); kFS_app(ii,ss) = Sol(5); kSF_app(ii,ss) = Sol(6); T2S_app(ii,ss) = Sol(7); T2F_app(ii,ss) = Sol(8);
        %end
        toc
        
    end
    
    T1S_final(ss) = mean(nonzeros(T1S_app(:,ss))); T1F_final(ss) = mean(nonzeros(T1F_app(:,ss))); T2S_final(ss) = mean(nonzeros(T2S_app(:,ss))); T2F_final(ss) = mean(nonzeros(T2F_app(:,ss))); M0F_final(ss) = mean(nonzeros(M0F_app(:,ss))); M0S_final(ss) = mean(nonzeros(M0S_app(:,ss))); kFS_final(ss) = mean(nonzeros(kFS_app(:,ss))); kSF_final(ss) = mean(nonzeros(kSF_app(:,ss)));
    
    % Calculate accuracy and precision of parameter estimates from fmincon - NOTE USE OF MWF.
    ErrorTP_T1S(ss) = abs(((T1S_final(ss) - T1_S)/T1_S)*100); ErrorTP_T2S(ss) = abs(((T2S_final(ss) - T2_S)/T2_S)*100);
    ErrorTP_T1F(ss) = abs(((T1F_final(ss) - T1_F)/T1_F)*100); ErrorTP_T2F(ss) = abs(((T2F_final(ss) - T2_F)/T2_F)*100);
    ErrorTP_M0F(ss) = abs(((M0F_final(ss) - M0_F)/M0_F)*100); ErrorTP_kFS(ss) = abs(((kFS_final(ss) - k_FS)/k_FS)*100);
    PrecisionTP_T1S(ss) = std(T1S_app(:,ss))/(T1S_final(ss)); PrecisionTP_T2S(ss) = std(T2S_app(:,ss))/(T2S_final(ss));
    PrecisionTP_T1F(ss) = std(T1F_app(:,ss))/(T1F_final(ss)); PrecisionTP_T2F(ss) = std(T2F_app(:,ss))/(T2F_final(ss));
    PrecisionTP_M0F(ss) = std(M0F_app(:,ss))/(M0F_final(ss)); PrecisionTP_kFS(ss) = std(kFS_app(:,ss))/mean(kFS_final(ss));
    
    %% Plot histograms of parameter estimates.
    
%     figure(2)
%     subplot(2,3,1); His3 = histogram(T1S_app(:,ss));
%     hold on; vline(T1_S,'k','GT')
%     xlabel('Apparent T_{1S}'); ylabel('Count');
%     His3.FaceColor = 'y'; His3.EdgeColor = 'k';
%     subplot(2,3,2); His4 = histogram(T1F_app(:,ss));
%     hold on; vline(T1_F,'k','GT')
%     xlabel('Apparent T_{1F}'); ylabel('Count');
%     His4.FaceColor = 'm'; His4.EdgeColor = 'k';
%     subplot(2,3,3); His5 = histogram(M0F_app(:,ss));
%     hold on; vline(M0_F(ss),'k','GT')
%     xlabel('Apparent M0_{F}'); ylabel('Count');
%     His5.FaceColor = 'c'; His5.EdgeColor = 'k';
%     subplot(2,3,4); His6 = histogram(kFS_app(:,ss));
%     hold on; vline(k_FS,'k','GT')
%     xlabel('Apparent k_{FS}'); ylabel('Count');
%     His6.FaceColor = 'r'; His6.EdgeColor = 'k';
%     subplot(2,3,5); His7 = histogram(T2S_app(:,ss));
%     hold on; vline(T2_S,'k','GT')
%     xlabel('Apparent T_{2S}'); ylabel('Count');
%     His7.FaceColor = 'g'; His7.EdgeColor = 'k';
%     subplot(2,3,6); His8 = histogram(T2F_app(:,ss));
%     hold on; vline(T2_F,'k','GT')
%     xlabel('Apparent T_{2F}'); ylabel('Count');
%     His8.FaceColor = 'b'; His8.EdgeColor = 'k';
        
    %% Signal Plots. Change SRC and two-stage lines for each sequence and functions if 3-pool fit used.
    
    % Calulcate parameter predictions from derived expressions.
    R1_F = 1/T1_F; R1_B = 1/T1_B; R1_S = 1/T1_S; dMzB = 0;
    G = 1.4e-5; Gamma = 2 * pi * 42.57747892e6; B1 = 13e-6;
    T_RF = deg2rad(50)/(Gamma * B1); W = (pi/T_RF) * (Gamma * B1)^2 * T_RF * G * (T_RF/TR_SPGR);
    % From symbolic maths.
    kSF_Prediction(ss) = (M0_F*k_FS)/M0_S + (M0_F*k_FB*k_SB)/(M0_B*W + M0_F*k_FB + M0_S*k_SB + M0_B*R1_B);
    kFS_Prediction(ss) = k_FS + (M0_S*k_FB*k_SB)/(M0_B*W + M0_F*k_FB + M0_S*k_SB + M0_B*R1_B);
    R1F_Prediction(ss) = (M0_B*R1_B*R1_F + M0_B*R1_F*W + M0_B*R1_B*k_FB + M0_F*R1_F*k_FB + M0_S*R1_F*k_SB + M0_B*W*k_FB)/(M0_B*W + M0_F*k_FB + M0_S*k_SB + M0_B*R1_B);
    R1S_Prediction(ss) = (M0_B*R1_B*R1_S + M0_B*R1_S*W + M0_F*R1_S*k_FB + M0_B*R1_B*k_SB + M0_S*R1_S*k_SB + M0_B*W*k_SB)/(M0_B*W + M0_F*k_FB + M0_S*k_SB + M0_B*R1_B);
    M0F_Prediction(ss) = (M0_F*(M0_B*R1_B*R1_F - dMzB*k_FB + M0_B*R1_F*W + M0_B*R1_B*k_FB + M0_F*R1_F*k_FB + M0_S*R1_F*k_SB))/(M0_B*R1_B*R1_F + M0_B*R1_F*W + M0_B*R1_B*k_FB + M0_F*R1_F*k_FB + M0_S*R1_F*k_SB + M0_B*W*k_FB);
    M0S_Prediction(ss) = (M0_S*(M0_B*R1_B*R1_S - dMzB*k_SB + M0_B*R1_S*W + M0_F*R1_S*k_FB + M0_B*R1_B*k_SB + M0_S*R1_S*k_SB))/(M0_B*R1_B*R1_S + M0_B*R1_S*W + M0_F*R1_S*k_FB + M0_B*R1_B*k_SB + M0_S*R1_S*k_SB + M0_B*W*k_SB);
    % From algebra.
    k_BS = k_SB*M0_S/M0_B; k_BF = k_FB*M0_F/M0_B; k_SF = k_FS*M0_F/M0_S;
    kFS_Prediction2(ss) = k_FS + ((k_BS*k_FB)/(R1_B + k_BS + k_BF + W));
    kSF_Prediction2(ss) = k_SF + ((k_BF*k_SB)/(R1_B + k_BS + k_BF + W));
    R1F_Prediction2(ss) = R1_F + k_FB - ((k_BS*k_FB + k_BF*k_FB)/(R1_B + k_BS + k_BF + W));    
    R1S_Prediction2(ss) = R1_S + k_SB - ((k_BF*k_SB + k_BS*k_SB)/(R1_B + k_BS + k_BF + W));
    M0F_Prediction2(ss) = (M0_F*R1_F + (((M0_B*k_BF*R1_B)-dMzB*k_BF)/(R1_B + k_BS + k_BF + W)))/R1F_Prediction2(ss);
    M0S_Prediction2(ss) = (M0_S*R1_S + (((M0_B*k_BS*R1_B)-dMzB*k_BS)/(R1_B + k_BS + k_BF + W)))/R1S_Prediction2(ss);
    
    figure(1)
    subplot(3,1,1)
    plot(rad2deg(FA_SPGR), (SPGR_steady_state_M0(FA_SPGR,TR_SPGR,'T1_S',T1S_final(ss),'T1_F',T1F_final(ss),'M0_F',M0F_final(ss),'M0_S',M0S_final(ss),'k_FS',kFS_final(ss),'k_SF',kSF_final(ss))),'k--','LineWidth',2)
    hold on
    plot(rad2deg(FA_SPGR), SPGR_Data, 'bo','LineWidth',2)
    %plot(rad2deg(FA_SPGR), abs(SPGR_Data_Noisy), 'bo','LineWidth',2)
    plot(rad2deg(FA_SPGR),SPGR_steady_state_M0(FA_SPGR,TR_SPGR,'T1_S',(1/R1S_Prediction(ss)),'T1_F',(1/R1F_Prediction(ss)),'M0_F',M0F_Prediction(ss),'M0_S',M0S_Prediction(ss),'k_FS',kFS_Prediction(ss),'k_SF',kSF_Prediction(ss)), 'r--','LineWidth', 2)
    plot(rad2deg(FA_SPGR),SPGR_steady_state_M0(FA_SPGR,TR_SPGR,'T1_S',(1/R1S_Prediction2(ss)),'T1_F',(1/R1F_Prediction2(ss)),'M0_F',M0F_Prediction2(ss),'M0_S',M0S_Prediction2(ss),'k_FS',kFS_Prediction2(ss),'k_SF',kSF_Prediction2(ss)), 'mo','LineWidth', 2)
    ll = legend('Simulation', 'GT', 'Prediction', 'Prediction2');
    ll.FontSize = 14;
    xlabel('FA [deg]','FontSize',14); ylabel('SPGR Signal','Fontsize',14)
    subplot(3,1,2)
    plot(rad2deg(FA_SSFP), (SSFP_steady_state_M0(FA_SSFP,TR_SSFP,'T1_S',T1S_final(ss),'T2_S',T2S_final(ss),'T1_F',T1F_final(ss),'T2_F',T2F_final(ss),'M0_F',M0F_final(ss),'M0_S',M0S_final(ss),'k_FS',kFS_final(ss),'k_SF',kSF_final(ss))),'k--','LineWidth',2)
    hold on
    plot(rad2deg(FA_SSFP), SSFP_Data,'bo','LineWidth',2)
    %plot(rad2deg(FA_SSFP), abs(SSFP_Data_Noisy),'bo','LineWidth',2)
    plot(rad2deg(FA_SSFP), SSFP_steady_state_M0(FA_SSFP, TR_SSFP, 'T1_S',(1/R1S_Prediction(ss)),'T2_S',T2_S,'T1_f',(1/R1F_Prediction(ss)),'T2_F',T2_F,'M0_F',M0F_Prediction(ss),'M0_S',M0S_Prediction(ss),'k_FS',kFS_Prediction(ss),'k_SF',kSF_Prediction(ss)), 'r--','LineWidth',2)
    plot(rad2deg(FA_SSFP), SSFP_steady_state_M0(FA_SSFP, TR_SSFP, 'T1_S',(1/R1S_Prediction2(ss)),'T2_S',T2_S,'T1_f',(1/R1F_Prediction2(ss)),'T2_F',T2_F,'M0_F',M0F_Prediction2(ss),'M0_S',M0S_Prediction2(ss),'k_FS',kFS_Prediction2(ss),'k_SF',kSF_Prediction2(ss)), 'mo','LineWidth',2)
    ll = legend('Simulation', 'GT', 'Prediction', 'Prediction2');
    ll.FontSize = 14;
    xlabel('FA [deg]','FontSize',14); ylabel('bSSFP_{0} Signal','Fontsize',14)
    subplot(3,1,3)
    plot(rad2deg(FA_SSFP), (SSFP_steady_state_180_M0(FA_SSFP,TR_SSFP,'T1_S',T1S_final(ss),'T2_S',T2S_final(ss),'T1_F',T1F_final(ss),'T2_F',T2F_final(ss),'M0_F',M0F_final(ss),'M0_S',M0S_final(ss),'k_FS',kFS_final(ss),'k_SF',kSF_final(ss))),'k--','LineWidth',2)
    hold on
    plot(rad2deg(FA_SSFP), SSFP_Data_180,'bo','LineWidth',2)
    %plot(rad2deg(FA_SSFP), abs(SSFP_Data_Noisy_180),'bo','LineWidth',2)
    plot(rad2deg(FA_SSFP), SSFP_steady_state_180_M0(FA_SSFP, TR_SSFP, 'T1_S',(1/R1S_Prediction(ss)),'T2_S',T2_S,'T1_F',(1/R1F_Prediction(ss)),'T2_F',T2_F,'M0_F',M0F_Prediction(ss),'M0_S',M0S_Prediction(ss),'k_FS',kFS_Prediction(ss),'k_SF',kSF_Prediction(ss)), 'r--','LineWidth',2)
    plot(rad2deg(FA_SSFP), SSFP_steady_state_180_M0(FA_SSFP, TR_SSFP, 'T1_S',(1/R1S_Prediction2(ss)),'T2_S',T2_S,'T1_F',(1/R1F_Prediction2(ss)),'T2_F',T2_F,'M0_F',M0F_Prediction2(ss),'M0_S',M0S_Prediction2(ss),'k_FS',kFS_Prediction2(ss),'k_SF',kSF_Prediction2(ss)), 'mo','LineWidth',2)
    ll = legend('Simulation', 'GT', 'Prediction', 'Prediction2');
    ll.FontSize = 14;
    xlabel('FA [deg]','FontSize',14); ylabel('bSSFP_{0} Signal','FontSize',14)
    xlabel('FA [deg]','FontSize',14); ylabel('bSSFP_{180} Signal','FontSize',14)
    
end

%% Predictions.

% k_bs_New = zeros(length(f_b), 1); k_bf_New = zeros(length(f_b), 1); f_s_New = zeros(length(f_b), 1); ff_Prediction = zeros(length(f_b), 1);
% 
% for tt = 1:length(f_b)
%     
%     f_s_New(tt) = 1 - (f_f + f_b(tt)); k_bs_New(tt) = (f_s_New(tt) * k_sb)/f_b(tt) ; k_bf_New(tt) = (f_f * k_fb)/f_b(tt);
%     R1_f = 1/T1_f; R1_b = 1/T1_b;
%     G = 1.4e-5; Gamma = 2 * pi * 42.57747892e6; B1 = 13e-6;
%     T_RF = deg2rad(50)/(Gamma * B1); W = (pi/T_RF) * (Gamma * B1)^2 * T_RF * G;
%     ff_Prediction(tt) = (R1_f * f_f + ((R1_b * f_b(tt) * k_bf_New(tt))/(R1_b + k_bf_New(tt) + k_bs_New(tt) + W)))/(R1_f - k_fb * (((k_bf_New(tt) + k_bs_New(tt))/(R1_b + k_bf_New(tt) + k_bs_New(tt) + W))-1));
%     
% end
% 
% figure(100)
% subplot(1,3,1)
% x = f_b; y = ff_final; errors = PrecisionTP_ff;
% errorbarxy(x, y, zeros(length(PrecisionTP_ff)), errors, {'.','k','k'});
% hold on
% errorbar(x, y, errors, 'o','Color','k','MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r')
% H = lsline;
% set(H, 'color', 'k', 'LineStyle','--')
% xlim([0 0.3]); ylim([0.15 0.4])
% xlabel('Bound-Pool Fraction','FontSize',14); ylabel('Simulated MWF','FontSize',14)
% 
% subplot(1,3,2)
% plot(f_b, ff_Prediction, 'ko', 'LineWidth', 1.5)
% xlabel('Bound-Pool Fraction','FontSize',14); ylabel('Predicted MWF','FontSize',14)
% 
% subplot(1,3,3)
% plot(ff_final, ff_Prediction, 'ro', 'LineWidth', 1.5)
% xlabel('Simulated MWF', 'FontSize', 14); ylabel('Predicted MWF', 'FontSize', 14)