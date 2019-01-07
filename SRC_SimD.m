%%% Runs SRC WITH exchange and B2. %%%

function [Sol, Bounds, Subsets] = SRC_SimD(Trials, Iterations, N, Runs, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data)
Bounds = zeros(Iterations,12); Subsets = zeros(5,6);

% Initialise loop variables.
T1S_app = zeros(Runs,1); M0F_app = zeros(Runs,1); kFS_app = zeros(Runs,1); T2S_app = zeros(Runs,1); T1F_app = zeros(Runs,1); T2F_app = zeros(Runs,1);

Candidates = zeros((length(FA_SPGR) + length(FA_SSFP0) + length(FA_SSFP180)),Trials);
Mss_SPGR = zeros(length(FA_SPGR),Trials); Mss_SSFP180 = zeros(length(FA_SSFP180),Trials); Mss_SSFP = zeros(length(FA_SSFP0),Trials);
SSFP_Signals = zeros((length(FA_SSFP0)+length(FA_SSFP180)),Trials); SRSQ = zeros(Trials,1);

T1S_Sub = zeros(N,1); T1F_Sub = zeros(N,1); T2S_Sub = zeros(N,1); T2F_Sub = zeros(N,1); M0F_Sub = zeros(N,1); kFS_Sub = zeros(N,1);
T1S_SRC = zeros(5,1); T1F_SRC = zeros(5,1); T2S_SRC = zeros(5,1); T2F_SRC = zeros(5,1); M0F_SRC = zeros(5,1); kFS_SRC = zeros(5,1);

T1S_LB = zeros(Iterations,1); T1S_UB = zeros(Iterations,1);
T2S_LB = zeros(Iterations,1); T2S_UB = zeros(Iterations,1);
T1F_LB = zeros(Iterations,1); T1F_UB = zeros(Iterations,1);
T2F_LB = zeros(Iterations,1); T2F_UB = zeros(Iterations,1);
M0F_LB = zeros(Iterations,1); M0F_UB = zeros(Iterations,1);
kFS_LB = zeros(Iterations,1); kFS_UB = zeros(Iterations,1);

for nn = 1:Runs
    
    % Bound set B2.
    T1F_LB(1) = 0.3; T1F_UB(1) = 0.65; T1S_LB(1) = 0.9; T1S_UB(1) = 5;
    T2F_LB(1) = 0.001; T2F_UB(1) = 0.03; T2S_LB(1) = 0.05; T2S_UB(1) = 0.165;
    M0F_LB(1) = 0; M0F_UB(1) = 0.35; kFS_LB(1) = (1/0.6); kFS_UB(1) = 40;
    
    for it = 1:Iterations
                
        if it == 1
            
            % Randomly choose parameter values from uniform distribution.
            T1S_Pick = (T1S_UB(it) - T1S_LB(it)) .* rand(Trials,1) + T1S_LB(it); T2S_Pick = (T2S_UB(it) - T2S_LB(it)) .* rand(Trials,1) + T2S_LB(it);
            T1F_Pick = (T1F_UB(it) - T1F_LB(it)) .* rand(Trials,1) + T1F_LB(it); T2F_Pick = (T2F_UB(it) - T2F_LB(it)) .* rand(Trials,1) + T2F_LB(it);
            M0F_Pick = (M0F_UB(it) - M0F_LB(it)) .* rand(Trials,1) + M0F_LB(it); kFS_Pick = (kFS_UB(it) - kFS_LB(it)) .* rand(Trials,1) + kFS_LB(it);
            
        else
            
            % Randomly choose parameter values from normal distribution.
            T1F_Pick = mean(T1F_Sub) + (std(T1F_Sub) .* trandn((((T1F_LB(it)*ones(Trials,1))-mean(T1F_Sub))./std(T1F_Sub)), (((T1F_UB(it)*ones(Trials,1))-mean(T1F_Sub))./std(T1F_Sub))));
            T1S_Pick = mean(T1S_Sub) + (std(T1S_Sub) .* trandn((((T1S_LB(it)*ones(Trials,1))-mean(T1S_Sub))./std(T1S_Sub)), (((T1S_UB(it)*ones(Trials,1))-mean(T1S_Sub))./std(T1S_Sub))));
            T2F_Pick = mean(T2F_Sub) + (std(T2F_Sub) .* trandn((((T2F_LB(it)*ones(Trials,1))-mean(T2F_Sub))./std(T2F_Sub)), (((T2F_UB(it)*ones(Trials,1))-mean(T2F_Sub))./std(T2F_Sub))));
            T2S_Pick = mean(T2S_Sub) + (std(T2S_Sub) .* trandn((((T2S_LB(it)*ones(Trials,1))-mean(T2S_Sub))./std(T2S_Sub)), (((T2S_UB(it)*ones(Trials,1))-mean(T2S_Sub))./std(T2S_Sub))));
            M0F_Pick = mean(M0F_Sub) + (std(M0F_Sub) .* trandn((((M0F_LB(it)*ones(Trials,1))-mean(M0F_Sub))./std(M0F_Sub)), (((M0F_UB(it)*ones(Trials,1))-mean(M0F_Sub))./std(M0F_Sub))));
            kFS_Pick = mean(kFS_Sub) + (std(kFS_Sub) .* trandn((((kFS_LB(it)*ones(Trials,1))-mean(kFS_Sub))./std(kFS_Sub)), (((kFS_UB(it)*ones(Trials,1))-mean(kFS_Sub))./std(kFS_Sub))));
            
        end
        
        for qq = 1:Trials
            
            [Mss_SPGR(:,qq), Mss_SSFP(:,qq), Mss_SSFP180(:,qq)] = MEX_mcDESPOT_B0_DE_DiffFAs(T1F_Pick(qq), T1S_Pick(qq), T2F_Pick(qq), T2S_Pick(qq), kFS_Pick(qq), M0F_Pick(qq), TR_SPGR, TR_SSFP, 0, rad2deg(FA_SPGR), rad2deg(FA_SSFP0), rad2deg(FA_SSFP180));
                        
            SSFP_Signals(:,qq) = [Mss_SSFP(:,qq); Mss_SSFP180(:,qq)];

            Candidates(:,qq) = [Mss_SPGR(:,qq)./mean(Mss_SPGR(:,qq)); SSFP_Signals(:,qq)./mean(SSFP_Signals(:,qq))];  
            SRSQ(qq) = norm(Candidates(:,qq) - Data).^2; 
            
        end
        
        % Sort residuals in ascending order to determine indices of N-lowest residual combinations.
        [~, Index] = sort(SRSQ);
        
        if it < Iterations
            
            % Locate indices in existing parameter sets.
            for pp = 1:N
                
                T1S_Sub(pp) = T1S_Pick(Index(pp)); T1F_Sub(pp) = T1F_Pick(Index(pp));
                T2S_Sub(pp) = T2S_Pick(Index(pp)); T2F_Sub(pp) = T2F_Pick(Index(pp));
                M0F_Sub(pp) = M0F_Pick(Index(pp)); kFS_Sub(pp) = kFS_Pick(Index(pp));
                
            end
            
            % Extract new search-space bounds from each subset.
            T1S_LB(it+1) = min(T1S_Sub) ; T1S_UB(it+1) = max(T1S_Sub); T1F_LB(it+1) = min(T1F_Sub) ; T1F_UB(it+1) = max(T1F_Sub);
            T2S_LB(it+1) = min(T2S_Sub) ; T2S_UB(it+1) = max(T2S_Sub); T2F_LB(it+1) = min(T2F_Sub) ; T2F_UB(it+1) = max(T2F_Sub);
            M0F_LB(it+1) = min(M0F_Sub) ; M0F_UB(it+1) = max(M0F_Sub); kFS_LB(it+1) = min(kFS_Sub) ; kFS_UB(it+1) = max(kFS_Sub);
            
            % Expand search-space slightly to avoid 'inadvertent over-contraction'.
            T1S_expn = (T1S_UB(it+1) - T1S_LB(it+1))/N; T1F_expn = (T1F_UB(it+1) - T1F_LB(it+1))/N;
            T2S_expn = (T2S_UB(it+1) - T2S_LB(it+1))/N; T2F_expn = (T2F_UB(it+1) - T2F_LB(it+1))/N;
            M0F_expn = (M0F_UB(it+1) - M0F_LB(it+1))/N; kFS_expn = (kFS_UB(it+1) - kFS_LB(it+1))/N;
            
            T1S_LB(it+1) = T1S_LB(it+1) - T1S_expn; T1S_UB(it+1) = T1S_UB(it+1) + T1S_expn;
            T1F_LB(it+1) = T1F_LB(it+1) - T1F_expn; T1F_UB(it+1) = T1F_UB(it+1) + T1F_expn;
            T2S_LB(it+1) = T2S_LB(it+1) - T2S_expn; T2S_UB(it+1) = T2S_UB(it+1) + T2S_expn;
            T2F_LB(it+1) = T2F_LB(it+1) - T2F_expn; T2F_UB(it+1) = T2F_UB(it+1) + T2F_expn;
            M0F_LB(it+1) = M0F_LB(it+1) - M0F_expn; M0F_UB(it+1) = M0F_UB(it+1) + M0F_expn;
            kFS_LB(it+1) = kFS_LB(it+1) - kFS_expn; kFS_UB(it+1) = kFS_UB(it+1) + kFS_expn;
            
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
            
            if (((max(T1S_Sub)-min(T1S_Sub))/max(T1S_Sub)) < 0.01 && ((max(T1F_Sub)-min(T1F_Sub))/max(T1F_Sub)) < 0.01 && ((max(M0F_Sub)-min(M0F_Sub))/max(M0F_Sub)) < 0.01 && ((max(T2S_Sub)-min(T2S_Sub))/max(T2S_Sub)) < 0.01 && ((max(T2F_Sub)-min(T2F_Sub))/max(T2F_Sub)) < 0.01 && ((max(kFS_Sub)-min(kFS_Sub))/max(kFS_Sub)) < 0.01)
               break
            end
            
        end
        
    end
    
    % Return result as mean of top-five solutions.
    for pp = 1:5
        T1S_SRC(pp) = T1S_Pick(Index(pp)); T1F_SRC(pp) = T1F_Pick(Index(pp));
        T2S_SRC(pp) = T2S_Pick(Index(pp)); T2F_SRC(pp) = T2F_Pick(Index(pp));
        M0F_SRC(pp) = M0F_Pick(Index(pp)); kFS_SRC(pp) = kFS_Pick(Index(pp));
    end
        
    T1S_app(nn) = mean(T1S_SRC); T1F_app(nn) = mean(T1F_SRC);
    T2S_app(nn) = mean(T2S_SRC); T2F_app(nn) = mean(T2F_SRC);
    M0F_app(nn) = mean(M0F_SRC); kFS_app(nn) = mean(kFS_SRC);
    
end

Sol(:,:) = [T1S_app(:), T1F_app(:), T2S_app(:), T2F_app(:), M0F_app(:), kFS_app(:)];
Bounds(:,:) = [T1S_LB(:), T1S_UB(:), T1F_LB(:), T1F_UB(:), T2S_LB(:), T2S_UB(:), T2F_LB(:), T2F_UB(:), M0F_LB(:), M0F_UB(:), kFS_LB(:), kFS_UB(:)];
Subsets(:,:) = [T1S_SRC, T1F_SRC, T2S_SRC, T2F_SRC, M0F_SRC, kFS_SRC];

end
