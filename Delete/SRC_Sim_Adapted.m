% SRC for pre-corrected FAs and pre-concatenated signals.

function [Sol,Bounds] = SRC_Sim_Adapted(Trials, Iterations, N, Runs, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data) %FA_SSFP0, FA_SSFP90, FA_SSFP270, 

rng('shuffle');

% Initialise loop variables.
T1S_Final = zeros(Runs,1); T1F_Final = zeros(Runs,1); 
T2S_Final = zeros(Runs,1); T2F_Final = zeros(Runs,1);
M0F_Final = zeros(Runs,1); Tau_Final = zeros(Runs,1);
%Delta_SRC = zeros(Runs,1);

Candidates = zeros((length(FA_SPGR) + length(FA_SSFP0) + length(FA_SSFP180)),Trials); %length(FA_SSFP90) + length(FA_SSFP270)
Mss_SPGR = zeros(length(FA_SPGR),Trials); Mss_SSFP180 = zeros(length(FA_SSFP180),Trials); Mss_SSFP = zeros(length(FA_SSFP0),Trials); %Mss_SSFP90 = zeros(length(FA_SSFP90),Trials); Mss_SSFP270 = zeros(length(FA_SSFP270),Trials);
SRSQ = zeros(Trials,1);
T1S_Sub = zeros(N,1); T1F_Sub = zeros(N,1); 
T2S_Sub = zeros(N,1); T2F_Sub = zeros(N,1); 
M0F_Sub = zeros(N,1); Tau_Sub = zeros(N,1);
%Delta_Sub = zeros(N,1);

T1S_LB = zeros(Iterations,1); T1S_UB = zeros(Iterations,1);
T2S_LB = zeros(Iterations,1); T2S_UB = zeros(Iterations,1);
T1F_LB = zeros(Iterations,1); T1F_UB = zeros(Iterations,1);
T2F_LB = zeros(Iterations,1); T2F_UB = zeros(Iterations,1);
M0F_LB = zeros(Iterations,1); M0F_UB = zeros(Iterations,1);
Tau_LB = zeros(Iterations,1); Tau_UB = zeros(Iterations,1);
%Delta_LB = zeros(Iterations+1,1); Delta_UB = zeros(Iterations+1,1);

for nn = 1:Runs
    % Hurley's Bounds.
    T1S_LB(1) = 0.5; T1S_UB(1) = 1.5; T2S_LB(1) = 0.05; T2S_UB(1) = 0.165;
    T1F_LB(1) = 0.3; T1F_UB(1) = 0.65; T2F_LB(1) = 0.001; T2F_UB(1) = 0.03;
    M0F_LB(1) = 0; M0F_UB(1) = 0.35; Tau_LB(1) = 0.025; Tau_UB(1) = 0.6;
    %Delta_LB(1) = -pi; Delta_UB(1) = pi;
    
    tic
    for it = 1:Iterations
        disp(['Iteration ', num2str(it)])
        
        if it == 1
            
            % Randomly choose parameter values from uniform distribution.
            T1S_Pick = (T1S_UB(it) - T1S_LB(it)) .* rand([Trials 1]) + T1S_LB(it); T2S_Pick = (T2S_UB(it) - T2S_LB(it)) .* rand([Trials 1]) + T2S_LB(it);
            T1F_Pick = (T1F_UB(it) - T1F_LB(it)) .* rand([Trials 1]) + T1F_LB(it); T2F_Pick = (T2F_UB(it) - T2F_LB(it)) .* rand([Trials 1]) + T2F_LB(it);
            M0F_Pick = (M0F_UB(it) - M0F_LB(it)) .* rand([Trials 1]) + M0F_LB(it); Tau_Pick = (Tau_UB(it) - Tau_LB(it)) .* rand([Trials 1]) + Tau_LB(it);
            %Delta_Pick = (Delta_UB(it) - Delta_LB(it)) .* rand(Trials,1) + Delta_LB(it);
            
            %T1S_Pick(1) = (0.965+0.465)/2; T1F_Pick(1) = (0.965+0.465)/2;
            %T2S_Pick(1) = (0.090+0.012)/2; T2F_Pick(1) = (0.090+0.012)/2;
            
        else
            
            % Randomly choose parameter values from normal distribution.
            %T1F_Pick = mean(T1F_Sub) + (std(T1F_Sub) .* trandn((((T1F_LB(it)*ones(Trials,1))-mean(T1F_Sub))./std(T1F_Sub)), (((T1F_UB(it)*ones(Trials,1))-mean(T1F_Sub))./std(T1F_Sub))));
            %T1S_Pick = mean(T1S_Sub) + (std(T1S_Sub) .* trandn((((T1S_LB(it)*ones(Trials,1))-mean(T1S_Sub))./std(T1S_Sub)), (((T1S_UB(it)*ones(Trials,1))-mean(T1S_Sub))./std(T1S_Sub))));
            %T2F_Pick = mean(T2F_Sub) + (std(T2F_Sub) .* trandn((((T2F_LB(it)*ones(Trials,1))-mean(T2F_Sub))./std(T2F_Sub)), (((T2F_UB(it)*ones(Trials,1))-mean(T2F_Sub))./std(T2F_Sub))));
            %T2S_Pick = mean(T2S_Sub) + (std(T2S_Sub) .* trandn((((T2S_LB(it)*ones(Trials,1))-mean(T2S_Sub))./std(T2S_Sub)), (((T2S_UB(it)*ones(Trials,1))-mean(T2S_Sub))./std(T2S_Sub))));
            %M0F_Pick = mean(M0F_Sub) + (std(M0F_Sub) .* trandn((((M0F_LB(it)*ones(Trials,1))-mean(M0F_Sub))./std(M0F_Sub)), (((M0F_UB(it)*ones(Trials,1))-mean(M0F_Sub))./std(M0F_Sub))));
            %kFS_Pick = mean(kFS_Sub) + (std(kFS_Sub) .* trandn((((kFS_LB(it)*ones(Trials,1))-mean(kFS_Sub))./std(kFS_Sub)), (((kFS_UB(it)*ones(Trials,1))-mean(kFS_Sub))./std(kFS_Sub))));
            %Delta_Pick = mean(Delta_Sub) + (std(Delta_Sub) .* trandn((((Delta_LB(it)*ones(Trials,1))-mean(Delta_Sub))./std(Delta_Sub)), (((Delta_UB(it)*ones(Trials,1))-mean(Delta_Sub))./std(Delta_Sub))));
            
            T1F_tmp = mean(T1F_Sub) * ones([Trials 1]) + randn([Trials 1]) * std(T1F_Sub);
            T1S_tmp = mean(T1S_Sub) * ones([Trials 1]) + randn([Trials 1]) * std(T1S_Sub);
            T2F_tmp = mean(T2F_Sub) * ones([Trials 1]) + randn([Trials 1]) * std(T2F_Sub);
            T2S_tmp = mean(T2S_Sub) * ones([Trials 1]) + randn([Trials 1]) * std(T2S_Sub);
            Tau_tmp = mean(Tau_Sub) * ones([Trials 1]) + randn([Trials 1]) * std(Tau_Sub);
            M0F_tmp = mean(M0F_Sub) * ones([Trials 1]) + randn([Trials 1]) * std(M0F_Sub);
            
            T1F_size = size(T1F_tmp((T1F_tmp < T1F_LB(it))));
            T1F_tmp((T1F_tmp < T1F_LB(it))) = mean(T1F_Sub) * ones([T1F_size 1]) + randn([T1F_size 1]) * std(T1F_Sub);
            T1F_tmp((T1F_tmp < T1F_LB(it))) = T1F_LB(it);
            
            T1S_size = size(T1S_tmp((T1S_tmp < T1S_LB(it))));
            T1S_tmp((T1S_tmp < T1S_LB(it))) = mean(T1S_Sub) * ones([T1S_size 1]) + randn([T1S_size 1]) * std(T1S_Sub);
            T1S_tmp((T1S_tmp < T1S_LB(it))) = T1S_LB(it);
               
            T2F_size = size(T2F_tmp((T2F_tmp < T2F_LB(it))));
            T2F_tmp((T2F_tmp < T2F_LB(it))) = mean(T2F_Sub) * ones([T2F_size 1]) + randn([T2F_size 1]) * std(T2F_Sub);
            T2F_tmp((T2F_tmp < T2F_LB(it))) = T2F_LB(it);
            
            T2S_size = size(T2S_tmp((T2S_tmp < T2S_LB(it))));
            T2S_tmp((T2S_tmp < T2S_LB(it))) = mean(T2S_Sub) * ones([T2S_size 1]) + randn([T2S_size 1]) * std(T2S_Sub);
            T2S_tmp((T2S_tmp < T2S_LB(it))) = T2S_LB(it);
            
            Tau_size = size(Tau_tmp((Tau_tmp < Tau_LB(it))));
            Tau_tmp((Tau_tmp < Tau_LB(it))) = mean(Tau_Sub) * ones([Tau_size 1]) + randn([Tau_size 1]) * std(Tau_Sub);
            Tau_tmp((Tau_tmp < Tau_LB(it))) = Tau_LB(it);
            
            M0F_size = size(M0F_tmp((M0F_tmp < M0F_LB(it))));
            M0F_tmp((M0F_tmp < M0F_LB(it))) = mean(M0F_Sub) * ones([M0F_size 1]) + randn([M0F_size 1]) * std(M0F_Sub);
            M0F_tmp((M0F_tmp < M0F_LB(it))) = M0F_LB(it);
             
            T1F_size = size(T1F_tmp((T1F_tmp > T1F_UB(it))));
            T1F_tmp((T1F_tmp > T1F_UB(it))) = mean(T1F_Sub) * ones([T1F_size 1]) + randn([T1F_size 1]) * std(T1F_Sub);
            T1F_tmp((T1F_tmp > T1F_UB(it))) = T1F_UB(it);
 
            T1S_size = size(T1S_tmp((T1S_tmp > T1S_UB(it))));
            T1S_tmp((T1S_tmp > T1S_UB(it))) = mean(T1S_Sub) * ones([T1S_size 1]) + randn([T1S_size 1]) * std(T1S_Sub);
            T1S_tmp((T1S_tmp > T1S_UB(it))) = T1S_UB(it);

            T2F_size = size(T2F_tmp((T2F_tmp > T2F_UB(it))));
            T2F_tmp((T2F_tmp > T2F_UB(it))) = mean(T2F_Sub) * ones([T2F_size 1]) + randn([T2F_size 1]) * std(T2F_Sub);
            T2F_tmp((T2F_tmp > T2F_UB(it))) = T2F_UB(it);
 
            T2S_size = size(T2S_tmp((T2S_tmp > T2S_UB(it))));
            T2S_tmp((T2S_tmp > T2S_UB(it))) = mean(T2S_Sub) * ones([T2S_size 1]) + randn([T2S_size 1]) * std(T2S_Sub);
            T2S_tmp((T2S_tmp > T2S_UB(it))) = T2S_UB(it);
            
            Tau_size = size(Tau_tmp((Tau_tmp > Tau_UB(it))));
            Tau_tmp((Tau_tmp > Tau_UB(it))) = mean(Tau_Sub) * ones([Tau_size 1]) + randn([Tau_size 1]) * std(Tau_Sub);
            Tau_tmp((Tau_tmp > Tau_UB(it))) = Tau_UB(it);
 
            M0F_size = size(M0F_tmp((M0F_tmp > M0F_UB(it))));
            M0F_tmp((M0F_tmp > M0F_UB(it))) = mean(M0F_Sub) * ones([M0F_size 1]) + randn([M0F_size 1]) * std(M0F_Sub);
            M0F_tmp((M0F_tmp > M0F_UB(it))) = M0F_UB(it);            
            
            T1F_Pick = T1F_tmp; T1S_Pick = T1S_tmp; T2F_Pick = T2F_tmp; T2S_Pick = T2S_tmp; Tau_Pick = Tau_tmp; M0F_Pick = M0F_tmp;
            
        end
        
        for qq = 1:Trials
            
            % Generate a number of theoretical signals, normalise and concatenate.
            %Mss_SPGR(:,qq) = SPGR_SteadyState(FA_SPGR,TR_SPGR,'T1_S',T1S_Pick(qq),'T1_F',T1F_Pick(qq),'M0_F',M0F_Pick(qq),'k_FS',kFS_Pick(qq));
            %Mss_SSFP(:,qq) = SSFP_SteadyState(FA_SSFP0,TR_SSFP,0,'T1_S',T1S_Pick(qq),'T2_S',T2S_Pick(qq),'T1_F',T1F_Pick(qq),'T2_F',T2F_Pick(qq),'M0_F',M0F_Pick(qq),'k_FS',kFS_Pick(qq));
            %Mss_SSFP180(:,qq) = SSFP_SteadyState(FA_SSFP180,TR_SSFP,pi,'T1_S',T1S_Pick(qq),'T2_S',T2S_Pick(qq),'T1_F',T1F_Pick(qq),'T2_F',T2F_Pick(qq),'M0_F',M0F_Pick(qq),'k_FS',kFS_Pick(qq));          
            
            [Mss_SPGR(:,qq), Mss_SSFP(:,qq), Mss_SSFP180(:,qq)] = MEX_mcDESPOT_B0_DE(T1F_Pick(qq), T1S_Pick(qq), T2F_Pick(qq), T2S_Pick(qq), (1/Tau_Pick(qq)), M0F_Pick(qq), TR_SPGR, TR_SSFP, 0, rad2deg(FA_SPGR), rad2deg(FA_SSFP180));
            %[Mss_SPGR(:,qq), ~, Mss_SSFP180(:,qq)] = MEX_mcDESPOT_B0_DE(T1F_Pick(qq), T1S_Pick(qq), T2F_Pick(qq), T2S_Pick(qq), kFS_Pick(qq), M0F_Pick(qq), TR_SPGR, TR_SSFP, Delta_Pick(qq), rad2deg(FA_SPGR), rad2deg(FA_SSFP180));           
            %[~, Mss_SSFP90(:,qq), ~] = MEX_mcDESPOT_B0_DE(T1F_Pick(qq), T1S_Pick(qq), T2F_Pick(qq), T2S_Pick(qq), kFS_Pick(qq), M0F_Pick(qq), TR_SPGR, TR_SSFP, ((pi/2)+Delta_Pick(qq)), rad2deg(FA_SPGR), rad2deg(FA_SSFP90));
            %[~, Mss_SSFP270(:,qq), ~] = MEX_mcDESPOT_B0_DE(T1F_Pick(qq), T1S_Pick(qq), T2F_Pick(qq), T2S_Pick(qq), kFS_Pick(qq), M0F_Pick(qq), TR_SPGR, TR_SSFP, ((3*pi/2)+Delta_Pick(qq)), rad2deg(FA_SPGR), rad2deg(FA_SSFP270));
            
            %Candidates(:,qq) = [Mss_SPGR(:,qq)./mean(Mss_SPGR(:,qq)); Mss_SSFP(:,qq)./mean(Mss_SSFP(:,qq)); Mss_SSFP180(:,qq)./mean(Mss_SSFP180(:,qq))];  % ; Mss_SSFP90(:,qq)./mean(Mss_SSFP90(:,qq)); Mss_SSFP270(:,qq)./mean(Mss_SSFP270(:,qq))  Candidates(:,qq) = [Mss_SPGR(:,qq); Mss_SSFP(:,qq); Mss_SSFP180(:,qq)];
            Candidates(:,qq) = [Mss_SPGR(:,qq); Mss_SSFP(:,qq); Mss_SSFP180(:,qq)];
            SRSQ(qq) = norm(Candidates(:,qq) - Data).^2;
            
        end
        
        % Sort residuals in ascending order to determine indices of N-lowest residual combinations.
        [~, Index] = sort(SRSQ);
        
        if it < Iterations
        
            % Locate indices in existing parameter sets.
            for pp = 1:N
                
                T1S_Sub(pp) = T1S_Pick(Index(pp)); T1F_Sub(pp) = T1F_Pick(Index(pp));
                T2S_Sub(pp) = T2S_Pick(Index(pp)); T2F_Sub(pp) = T2F_Pick(Index(pp));
                M0F_Sub(pp) = M0F_Pick(Index(pp)); Tau_Sub(pp) = Tau_Pick(Index(pp));
                %Delta_Sub(pp) = Delta_Pick(Index(pp));
                
            end
            
            % Extract new search-space bounds from each subset.
            T1S_LB(it+1) = min(T1S_Sub) ; T1S_UB(it+1) = max(T1S_Sub); T1F_LB(it+1) = min(T1F_Sub) ; T1F_UB(it+1) = max(T1F_Sub);
            T2S_LB(it+1) = min(T2S_Sub) ; T2S_UB(it+1) = max(T2S_Sub); T2F_LB(it+1) = min(T2F_Sub) ; T2F_UB(it+1) = max(T2F_Sub);
            M0F_LB(it+1) = min(M0F_Sub) ; M0F_UB(it+1) = max(M0F_Sub); Tau_LB(it+1) = min(Tau_Sub) ; Tau_UB(it+1) = max(Tau_Sub);
            %Delta_LB(it+1) = min(Delta_Sub) ; Delta_UB(it+1) = max(Delta_Sub);
            
            % Expand search-space slightly to avoid 'inadvertent over-contraction'.
            T1S_expn = (T1S_UB(it+1) - T1S_LB(it+1))/N; T1F_expn = (T1F_UB(it+1) - T1F_LB(it+1))/N;
            T2S_expn = (T2S_UB(it+1) - T2S_LB(it+1))/N; T2F_expn = (T2F_UB(it+1) - T2F_LB(it+1))/N;
            M0F_expn = (M0F_UB(it+1) - M0F_LB(it+1))/N; Tau_expn = (Tau_UB(it+1) - Tau_LB(it+1))/N;
            %Delta_expn = (Delta_UB(it+1) - Delta_LB(it+1))/N;
            
            
            T1S_LB(it+1) = T1S_LB(it+1) - T1S_expn; T1S_UB(it+1) = T1S_UB(it+1) + T1S_expn;
            T1F_LB(it+1) = T1F_LB(it+1) - T1F_expn; T1F_UB(it+1) = T1F_UB(it+1) + T1F_expn;
            T2S_LB(it+1) = T2S_LB(it+1) - T2S_expn; T2S_UB(it+1) = T2S_UB(it+1) + T2S_expn;
            T2F_LB(it+1) = T2F_LB(it+1) - T2F_expn; T2F_UB(it+1) = T2F_UB(it+1) + T2F_expn;
            M0F_LB(it+1) = M0F_LB(it+1) - M0F_expn; M0F_UB(it+1) = M0F_UB(it+1) + M0F_expn;
            Tau_LB(it+1) = Tau_LB(it+1) - Tau_expn; Tau_UB(it+1) = Tau_UB(it+1) + Tau_expn;
            %Delta_LB(it+1) = Delta_LB(it+1) - Delta_expn; Delta_UB(it+1) = Delta_UB(it+1) + Delta_expn;
            
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
            if Tau_LB(it+1) < 0
                Tau_LB(it+1) = 0;
            end
            %if Delta_LB(it+1) < 0
            %   Delta_LB(it+1) = 0;
            %end
            %if (((max(T1S_Sub)-min(T1S_Sub))/max(T1S_Sub)) < 0.01 && ((max(T1F_Sub)-min(T1F_Sub))/max(T1F_Sub)) < 0.01 && ((max(kFS_Sub)-min(kFS_Sub))/max(kFS_Sub)) < 0.01 && ((max(M0F_Sub)-min(M0F_Sub))/max(M0F_Sub)) < 0.01 && ((max(T2S_Sub)-min(T2S_Sub))/max(T2S_Sub)) < 0.01 && ((max(Delta_Sub)-min(Delta_Sub))/max(Delta_Sub)) < 0.01 && ((max(T2F_Sub)-min(T2F_Sub))/max(T2F_Sub)) < 0.01)
            %   break
            %end
        end
        
    end
    toc
   
    for pp = 1:5
        
        T1S_SRC = T1S_Pick(Index(pp)); T1F_SRC = T1F_Pick(Index(pp));
        T2S_SRC = T2S_Pick(Index(pp)); T2F_SRC = T2F_Pick(Index(pp));
        M0F_SRC = M0F_Pick(Index(pp)); Tau_SRC = Tau_Pick(Index(pp));
        %Delta_SRC = Delta_Pick(Index(pp));
        
    end
    
    % UPDATE GUESS AGAIN BUT ONLY MEAN NEEDED. Return result as mean of top-five solutions.
    T1S_Final(nn) = mean(T1S_SRC); T1F_Final(nn) = mean(T1F_SRC); 
    T2S_Final(nn) = mean(T2S_SRC); T2F_Final(nn) = mean(T2F_SRC); 
    M0F_Final(nn) = mean(M0F_SRC); Tau_Final(nn) = mean(Tau_SRC); 
    %Delta_Final(nn) = mean(Delta_SRC);  
    Sol(:,:) = [T1S_Final(:), T1F_Final(:), T2S_Final(:), T2F_Final(:), M0F_Final(:), Tau_Final(:)]; %, Delta_Final(:)
    Bounds(:,:) = [T1F_LB(:), T1F_UB(:), T1S_LB(:), T1S_UB(:), T2F_LB(:), T2F_UB(:), T2S_LB(:), T2S_UB(:), M0F_LB(:), M0F_UB(:), Tau_LB(:), Tau_UB(:)];
    
end

end
