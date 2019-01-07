%%% SRC for in vivo data. Model INCLUDES exchange. %%%

function [T1S_final, T1F_final, T2S_final, T2F_final, kFS_final, M0F_final, Delta_final] = SRC_mcDESPOT_Wood(Coords, Trials, Iterations, N, FA_SPGR, FA_SSFP180, FA_SSFP0, TR_SPGR, TR_SSFP, Data)

T1S_final = zeros(length(Coords),1); T2S_final = zeros(length(Coords),1); M0F_final = zeros(length(Coords),1); kFS_final = zeros(length(Coords),1); T1F_final = zeros(length(Coords),1); T2F_final = zeros(length(Coords),1); Delta_final = zeros(length(Coords),1);

parfor mm = 1:length(Coords)
       
    % Initialise loop variables.
    Candidates = zeros((length(FA_SPGR(mm,:)) + length(FA_SSFP180(mm,:)) + length(FA_SSFP0(mm,:))),Trials);
    Mss_SPGR = zeros(length(FA_SPGR(mm,:)),Trials); Mss_SSFP0 = zeros(length(FA_SSFP0(mm,:)),Trials); Mss_SSFP180 = zeros(length(FA_SSFP180(mm,:)),Trials);
    SSFP_Signals = zeros((length(FA_SSFP0(mm,:))+length(FA_SSFP180(mm,:))),Trials);
    SRSQ = zeros(Trials,1);
    T1S_Sub = zeros(N,1); T2S_Sub = zeros(N,1); M0F_Sub = zeros(N,1); kFS_Sub = zeros(N,1); T1F_Sub = zeros(N,1); T2F_Sub = zeros(N,1); Delta_Sub = zeros(N,1);
    T1S_SRC = zeros(5,1); T2S_SRC = zeros(5,1); M0F_SRC = zeros(5,1); kFS_SRC = zeros(5,1); T1F_SRC = zeros(5,1); T2F_SRC = zeros(5,1); Delta_SRC = zeros(5,1);
    
    T1S_LB = zeros(Iterations,1); T1S_UB = zeros(Iterations,1);
    T2S_LB = zeros(Iterations,1); T2S_UB = zeros(Iterations,1);
    T1F_LB = zeros(Iterations,1); T1F_UB = zeros(Iterations,1);
    T2F_LB = zeros(Iterations,1); T2F_UB = zeros(Iterations,1);
    M0F_LB = zeros(Iterations,1); M0F_UB = zeros(Iterations,1);
    kFS_LB = zeros(Iterations,1); kFS_UB = zeros(Iterations,1);
    Delta_LB = zeros(Iterations,1); Delta_UB = zeros(Iterations,1);
    
    % Specify initial upper and lower bounds. B3
    T1S_LB(1) = 0.9; T1S_UB(1) = 1.5; T2S_LB(1) = 0.04; T2S_UB(1) = 0.15;
    T1F_LB(1) = 0.3; T1F_UB(1) = 0.8; T2F_LB(1) = 0.01; T2F_UB(1) = 0.03;
    M0F_LB(1) = 0.001; M0F_UB(1) = 0.35; kFS_LB(1) = (1/0.6); kFS_UB(1) = 40;
    Delta_LB(1) = 0; Delta_UB(1) = 2*pi;
    
    for it = 1:Iterations
        
        if it == 1
            
            % Randomly choose parameter values from uniform distribution.
            T1S_Pick = (T1S_UB(it) - T1S_LB(it)) .* rand(Trials,1) + T1S_LB(it); T1F_Pick = (T1F_UB(it) - T1F_LB(it)) .* rand(Trials,1) + T1F_LB(it);
            T2S_Pick = (T2S_UB(it) - T2S_LB(it)) .* rand(Trials,1) + T2S_LB(it); T2F_Pick = (T2F_UB(it) - T2F_LB(it)) .* rand(Trials,1) + T2F_LB(it);
            M0F_Pick = (M0F_UB(it) - M0F_LB(it)) .* rand(Trials,1) + M0F_LB(it); kFS_Pick = (kFS_UB(it) - kFS_LB(it)) .* rand(Trials,1) + kFS_LB(it);
            Delta_Pick = (Delta_UB(it) - Delta_LB(it)) .* rand(Trials,1) + Delta_LB(it);
            
        else
            
            T1F_Pick = mean(T1F_Sub) + (std(T1F_Sub) .* trandn((((T1F_LB(it)*ones(Trials,1))-mean(T1F_Sub))./std(T1F_Sub)), (((T1F_UB(it)*ones(Trials,1))-mean(T1F_Sub))./std(T1F_Sub))));
            T1S_Pick = mean(T1S_Sub) + (std(T1S_Sub) .* trandn((((T1S_LB(it)*ones(Trials,1))-mean(T1S_Sub))./std(T1S_Sub)), (((T1S_UB(it)*ones(Trials,1))-mean(T1S_Sub))./std(T1S_Sub))));
            T2F_Pick = mean(T2F_Sub) + (std(T2F_Sub) .* trandn((((T2F_LB(it)*ones(Trials,1))-mean(T2F_Sub))./std(T2F_Sub)), (((T2F_UB(it)*ones(Trials,1))-mean(T2F_Sub))./std(T2F_Sub))));
            T2S_Pick = mean(T2S_Sub) + (std(T2S_Sub) .* trandn((((T2S_LB(it)*ones(Trials,1))-mean(T2S_Sub))./std(T2S_Sub)), (((T2S_UB(it)*ones(Trials,1))-mean(T2S_Sub))./std(T2S_Sub))));
            M0F_Pick = mean(M0F_Sub) + (std(M0F_Sub) .* trandn((((M0F_LB(it)*ones(Trials,1))-mean(M0F_Sub))./std(M0F_Sub)), (((M0F_UB(it)*ones(Trials,1))-mean(M0F_Sub))./std(M0F_Sub))));
            kFS_Pick = mean(kFS_Sub) + (std(kFS_Sub) .* trandn((((kFS_LB(it)*ones(Trials,1))-mean(kFS_Sub))./std(kFS_Sub)), (((kFS_UB(it)*ones(Trials,1))-mean(kFS_Sub))./std(kFS_Sub))));
            Delta_Pick = mean(Delta_Sub) + (std(Delta_Sub) .* trandn((((Delta_LB(it)*ones(Trials,1))-mean(Delta_Sub))./std(Delta_Sub)), (((Delta_UB(it)*ones(Trials,1))-mean(Delta_Sub))./std(Delta_Sub))));
            
        end
        
        for qq = 1:Trials
            
            [Mss_SPGR(:,qq), Mss_SSFP0(:,qq), Mss_SSFP180(:,qq)] = MEX_mcDESPOT_B0_DE_DiffFAs(T1F_Pick(qq), T1S_Pick(qq), T2F_Pick(qq), T2S_Pick(qq), kFS_Pick(qq), M0F_Pick(qq), TR_SPGR, TR_SSFP, Delta_Pick(qq), FA_SPGR(mm,:), FA_SSFP0(mm,:), FA_SSFP180(mm,:));
            SSFP_Signals(:,qq) = [Mss_SSFP180(:,qq); Mss_SSFP0(:,qq)];
            Candidates(:,qq) = [Mss_SPGR(:,qq)./mean(Mss_SPGR(:,qq),1) ; SSFP_Signals(:,qq)./mean(SSFP_Signals(:,qq),1)];
            SRSQ(qq) = norm(Candidates(:,qq) - Data(:,mm))^2;
            
        end
        
        % Sort residuals in ascending order to determine indices of N-lowest residual combinations.
        [~, Index] = sort(SRSQ);
        
        if it < Iterations
            
            % Locate indices in existing parameter sets.
            for pp = 1:N
                
                T1S_Sub(pp) = T1S_Pick(Index(pp)); T1F_Sub(pp) = T1F_Pick(Index(pp));
                T2S_Sub(pp) = T2S_Pick(Index(pp)); T2F_Sub(pp) = T2F_Pick(Index(pp));
                M0F_Sub(pp) = M0F_Pick(Index(pp)); kFS_Sub(pp) = kFS_Pick(Index(pp));
                Delta_Sub(pp) = Delta_Pick(Index(pp));
                
            end
            
            % Extract new search-space bounds from each subset.
            T1S_LB(it+1) = min(T1S_Sub) ; T1S_UB(it+1) = max(T1S_Sub); T1F_LB(it+1) = min(T1F_Sub) ; T1F_UB(it+1) = max(T1F_Sub);
            M0F_LB(it+1) = min(M0F_Sub) ; M0F_UB(it+1) = max(M0F_Sub); kFS_LB(it+1) = min(kFS_Sub) ; kFS_UB(it+1) = max(kFS_Sub);
            T2S_LB(it+1) = min(T2S_Sub) ; T2S_UB(it+1) = max(T2S_Sub); T2F_LB(it+1) = min(T2F_Sub) ; T2F_UB(it+1) = max(T2F_Sub);
            Delta_LB(it+1) = min(Delta_Sub); Delta_UB(it+1) = max(Delta_Sub);
            
            % Expand search-space slightly to avoid 'inadvertent over-contraction'.
            T1S_expn = (T1S_UB(it+1) - T1S_LB(it+1))/N; T1F_expn = (T1F_UB(it+1) - T1F_LB(it+1))/N;
            T2S_expn = (T2S_UB(it+1) - T2S_LB(it+1))/N; T2F_expn = (T2F_UB(it+1) - T2F_LB(it+1))/N;
            M0F_expn = (M0F_UB(it+1) - M0F_LB(it+1))/N; kFS_expn = (kFS_UB(it+1) - kFS_LB(it+1))/N;
            Delta_expn = (Delta_UB(it+1) - Delta_LB(it+1))/N;
            
            T1S_LB(it+1) = T1S_LB(it+1) - T1S_expn; T1S_UB(it+1) = T1S_UB(it+1) + T1S_expn;
            T1F_LB(it+1) = T1F_LB(it+1) - T1F_expn; T1F_UB(it+1) = T1F_UB(it+1) + T1F_expn;
            T2S_LB(it+1) = T2S_LB(it+1) - T2S_expn; T2S_UB(it+1) = T2S_UB(it+1) + T2S_expn;
            T2F_LB(it+1) = T2F_LB(it+1) - T2F_expn; T2F_UB(it+1) = T2F_UB(it+1) + T2F_expn;
            M0F_LB(it+1) = M0F_LB(it+1) - M0F_expn; M0F_UB(it+1) = M0F_UB(it+1) + M0F_expn;
            kFS_LB(it+1) = kFS_LB(it+1) - kFS_expn; kFS_UB(it+1) = kFS_UB(it+1) + kFS_expn;
            Delta_LB(it+1) = Delta_LB(it+1) - Delta_expn; Delta_UB(it+1) = Delta_UB(it+1) + Delta_expn;
            
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
            if Delta_LB(it+1) < 0
                Delta_LB(it+1) = 0;
            end
            
            % Exit loop if difference between min and max of each parameter < 1%.
            if (((max(T1S_Sub)-min(T1S_Sub))/max(T1S_Sub)) < 0.01 && ((max(T1F_Sub)-min(T1F_Sub))/max(T1F_Sub)) < 0.01 && ((max(M0F_Sub)-min(M0F_Sub))/max(M0F_Sub)) < 0.01 && ((max(kFS_Sub)-min(kFS_Sub))/max(kFS_Sub)) < 0.01 && ((max(T2S_Sub)-min(T2S_Sub))/max(T2S_Sub)) < 0.01 && ((max(T2F_Sub)-min(T2F_Sub))/max(T2F_Sub)) < 0.01 && ((max(Delta_Sub)-min(Delta_Sub))/max(Delta_Sub)) < 0.01)
                break
            end
            
        end
        
    end
    
    % Return result as mean of top-five solutions.
    for pp = 1:5
        T1S_SRC(pp) = T1S_Pick(Index(pp)); T1F_SRC(pp) = T1F_Pick(Index(pp));
        T2S_SRC(pp) = T2S_Pick(Index(pp)); T2F_SRC(pp) = T2F_Pick(Index(pp));
        M0F_SRC(pp) = M0F_Pick(Index(pp)); kFS_SRC(pp) = kFS_Pick(Index(pp));
        Delta_SRC(pp) = Delta_Pick(Index(pp));
    end
    
    T1S_final(mm) = mean(T1S_SRC); T2S_final(mm) = mean(T2S_SRC); M0F_final(mm) = mean(M0F_SRC); kFS_final(mm) = mean(kFS_SRC); T1F_final(mm) = mean(T1F_SRC); T2F_final(mm) = mean(T2F_SRC); Delta_final(mm) = mean(Delta_SRC);
    
end

end
