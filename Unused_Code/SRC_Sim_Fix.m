% SRC for pre-corrected FAs and pre-concatenated signals.

function [Sol, Bounds, Subsets] = SRC_Sim_Fix(Trials, Iterations, N, Runs, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data) %FA_SSFP0, FA_SSFP90, FA_SSFP270,
Bounds = zeros(Iterations,14); Subsets = zeros(5,7);

% Initialise loop variables.
T1S_app = zeros(Runs,1); M0F_app = zeros(Runs,1); kFS_app = zeros(Runs,1); T2S_app = zeros(Runs,1); T1F_app = zeros(Runs,1); T2F_app = zeros(Runs,1); Delta_app = zeros(Runs,1);

Candidates = zeros((length(FA_SPGR) + length(FA_SSFP0) + length(FA_SSFP180)),Trials); %length(FA_SSFP90) + length(FA_SSFP270)
Mss_SPGR = zeros(length(FA_SPGR),Trials); Mss_SSFP180 = zeros(length(FA_SSFP180),Trials); Mss_SSFP = zeros(length(FA_SSFP0),Trials); %Mss_SSFP90 = zeros(length(FA_SSFP90),Trials); Mss_SSFP270 = zeros(length(FA_SSFP270),Trials);
SSFP_Signals = zeros((length(FA_SSFP0)+length(FA_SSFP180)),Trials); SRSQ = zeros(Trials,1);
T1S_Sub = zeros(N,1); T1F_Sub = zeros(N,1);
T2S_Sub = zeros(N,1); T2F_Sub = zeros(N,1);
M0F_Sub = zeros(N,1); kFS_Sub = zeros(N,1);
Delta_Sub = zeros(N,1);

T1S_LB = zeros(Iterations,1); T1S_UB = zeros(Iterations,1);
T2S_LB = zeros(Iterations,1); T2S_UB = zeros(Iterations,1);
T1F_LB = zeros(Iterations,1); T1F_UB = zeros(Iterations,1);
T2F_LB = zeros(Iterations,1); T2F_UB = zeros(Iterations,1);
M0F_LB = zeros(Iterations,1); M0F_UB = zeros(Iterations,1);
kFS_LB = zeros(Iterations,1); kFS_UB = zeros(Iterations,1);
Delta_LB = zeros(Iterations,1); Delta_UB = zeros(Iterations,1);

for nn = 1:Runs
    % Toby's Bounds.
    T1S_LB(1) = 0.9; T1S_UB(1) = 1.5; T2S_LB(1) = 0.04; T2S_UB(1) = 0.15;
    T1F_LB(1) = 0.3; T1F_UB(1) = 0.8; T2F_LB(1) = 0.01; T2F_UB(1) = 0.03;
    M0F_LB(1) = 0.001; M0F_UB(1) = 0.35; kFS_LB(1) = (1/0.6); kFS_UB(1) = 40;
    Delta_LB(1) = -pi; Delta_UB(1) = pi;
    % Zhang Bounds.
    %T1S_LB(1) = 0.7; T1S_UB(1) = 2.5; T2S_LB(1) = 0.075; T2S_UB(1) = 0.2;
    %T1F_LB(1) = 0.2; T1F_UB(1) = 0.5; T2F_LB(1) = 0.002; T2F_UB(1) = 0.045;
    %M0F_LB(1) = 1e-7; M0F_UB(1) = 0.30; kFS_LB(1) = 0.5; kFS_UB(1) = 20;
    %Delta_LB(1) = -pi; Delta_UB(1) = pi;
    % Bouhrara Bounds + Strict k Bounds.
    %T1S_LB(1) = 0.01; T1S_UB(1) = 3.5; T2S_LB(1) = 0.06; T2S_UB(1) = 0.2;
    %T1F_LB(1) = 0.01; T1F_UB(1) = 0.65; T2F_LB(1) = 0.001; T2F_UB(1) = 0.06;
    %M0F_LB(1) = 0; M0F_UB(1) = 0.45; kFS_LB(1) = 0.5; kFS_UB(1) = 20;
    %Delta_LB(1) = -pi; Delta_UB(1) = pi;
    % Deoni Bounds - older papers assume same bounds for F & S.
    %T1F_LB(1) = 0.3; T1F_UB(1) = 0.65; T1S_LB(1) = 0.9; T1S_UB(1) = 5;
    %T2F_LB(1) = 0.001; T2F_UB(1) = 0.03; T2S_LB(1) = 0.05; T2S_UB(1) = 0.165;
    %M0F_LB(1) = 0; M0F_UB(1) = 0.35; kFS_LB(1) = (1/0.6); kFS_UB(1) = 40;
    %Delta_LB(1) = 0; Delta_UB(1) = 2*pi;
    
    tic
    for it = 1:Iterations
        disp(['Iteration ', num2str(it)])
        
        if it == 1
            
            % Randomly choose parameter values from uniform distribution.
            T1S_Pick = (T1S_UB(it) - T1S_LB(it)) .* rand(Trials,1) + T1S_LB(it); T2S_Pick = (T2S_UB(it) - T2S_LB(it)) .* rand(Trials,1) + T2S_LB(it);
            T1F_Pick = (T1F_UB(it) - T1F_LB(it)) .* rand(Trials,1) + T1F_LB(it); T2F_Pick = (T2F_UB(it) - T2F_LB(it)) .* rand(Trials,1) + T2F_LB(it);
            M0F_Pick = (M0F_UB(it) - M0F_LB(it)) .* rand(Trials,1) + M0F_LB(it); kFS_Pick = (kFS_UB(it) - kFS_LB(it)) .* rand(Trials,1) + kFS_LB(it);
            Delta_Pick = (Delta_UB(it) - Delta_LB(it)) .* rand(Trials,1) + Delta_LB(it);
            
            %T1S_Pick(1) = (0.965+0.465)/2; T1F_Pick(1) = (0.965+0.465)/2;
            %T2S_Pick(1) = (0.090+0.012)/2; T2F_Pick(1) = (0.090+0.012)/2;
            
        else
            
            % Randomly choose parameter values from normal distribution.
            T1F_Pick = mean(T1F_Sub) + (std(T1F_Sub) .* trandn((((T1F_LB(it)*ones(Trials,1))-mean(T1F_Sub))./std(T1F_Sub)), (((T1F_UB(it)*ones(Trials,1))-mean(T1F_Sub))./std(T1F_Sub))));
            T1S_Pick = mean(T1S_Sub) + (std(T1S_Sub) .* trandn((((T1S_LB(it)*ones(Trials,1))-mean(T1S_Sub))./std(T1S_Sub)), (((T1S_UB(it)*ones(Trials,1))-mean(T1S_Sub))./std(T1S_Sub))));
            T2F_Pick = mean(T2F_Sub) + (std(T2F_Sub) .* trandn((((T2F_LB(it)*ones(Trials,1))-mean(T2F_Sub))./std(T2F_Sub)), (((T2F_UB(it)*ones(Trials,1))-mean(T2F_Sub))./std(T2F_Sub))));
            T2S_Pick = mean(T2S_Sub) + (std(T2S_Sub) .* trandn((((T2S_LB(it)*ones(Trials,1))-mean(T2S_Sub))./std(T2S_Sub)), (((T2S_UB(it)*ones(Trials,1))-mean(T2S_Sub))./std(T2S_Sub))));
            M0F_Pick = mean(M0F_Sub) + (std(M0F_Sub) .* trandn((((M0F_LB(it)*ones(Trials,1))-mean(M0F_Sub))./std(M0F_Sub)), (((M0F_UB(it)*ones(Trials,1))-mean(M0F_Sub))./std(M0F_Sub))));
            kFS_Pick = mean(kFS_Sub) + (std(kFS_Sub) .* trandn((((kFS_LB(it)*ones(Trials,1))-mean(kFS_Sub))./std(kFS_Sub)), (((kFS_UB(it)*ones(Trials,1))-mean(kFS_Sub))./std(kFS_Sub))));
            Delta_Pick = mean(Delta_Sub) + (std(Delta_Sub) .* trandn((((Delta_LB(it)*ones(Trials,1))-mean(Delta_Sub))./std(Delta_Sub)), (((Delta_UB(it)*ones(Trials,1))-mean(Delta_Sub))./std(Delta_Sub))));
            
        end
        
        for qq = 1:Trials
            
            [Mss_SPGR(:,qq), Mss_SSFP(:,qq), Mss_SSFP180(:,qq)] = MEX_mcDESPOT_B0_DE_DiffFAs(0.465, T1S_Pick(qq), 0.012, T2S_Pick(qq), 8, M0F_Pick(qq), TR_SPGR, TR_SSFP, Delta_Pick(qq), rad2deg(FA_SPGR), rad2deg(FA_SSFP0), rad2deg(FA_SSFP180));
            
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
                Delta_Sub(pp) = Delta_Pick(Index(pp));
                
            end
            
            % Extract new search-space bounds from each subset.
            T1S_LB(it+1) = min(T1S_Sub) ; T1S_UB(it+1) = max(T1S_Sub); T1F_LB(it+1) = min(T1F_Sub) ; T1F_UB(it+1) = max(T1F_Sub);
            T2S_LB(it+1) = min(T2S_Sub) ; T2S_UB(it+1) = max(T2S_Sub); T2F_LB(it+1) = min(T2F_Sub) ; T2F_UB(it+1) = max(T2F_Sub);
            M0F_LB(it+1) = min(M0F_Sub) ; M0F_UB(it+1) = max(M0F_Sub); kFS_LB(it+1) = min(kFS_Sub) ; kFS_UB(it+1) = max(kFS_Sub);
            Delta_LB(it+1) = min(Delta_Sub) ; Delta_UB(it+1) = max(Delta_Sub);
            
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
            %if Delta_LB(it+1) < 0; Delta_LB(it+1) = 0; end;
            
            %if (((max(T1S_Sub)-min(T1S_Sub))/max(T1S_Sub)) < 0.01 && ((max(T1F_Sub)-min(T1F_Sub))/max(T1F_Sub)) < 0.01 && ((max(kFS_Sub)-min(kFS_Sub))/max(kFS_Sub)) < 0.01 && ((max(M0F_Sub)-min(M0F_Sub))/max(M0F_Sub)) < 0.01 && ((max(T2S_Sub)-min(T2S_Sub))/max(T2S_Sub)) < 0.01 && ((max(Delta_Sub)-min(Delta_Sub))/max(Delta_Sub)) < 0.01 && ((max(T2F_Sub)-min(T2F_Sub))/max(T2F_Sub)) < 0.01)
            %   break
            %end
            
        end
        
    end
    toc
    % Return result as mean of top-five solutions.
    for pp = 1:5
        T1S_SRC = T1S_Pick(Index(pp)); T1F_SRC = T1F_Pick(Index(pp));
        T2S_SRC = T2S_Pick(Index(pp)); T2F_SRC = T2F_Pick(Index(pp));
        M0F_SRC = M0F_Pick(Index(pp)); kFS_SRC = kFS_Pick(Index(pp));
        Delta_SRC = Delta_Pick(Index(pp));
    end
        
    T1S_app(nn) = mean(T1S_SRC); T1F_app(nn) = mean(T1F_SRC);
    T2S_app(nn) = mean(T2S_SRC); T2F_app(nn) = mean(T2F_SRC);
    M0F_app(nn) = mean(M0F_SRC); kFS_app(nn) = mean(kFS_SRC);
    Delta_app(nn) = mean(Delta_SRC);
    
end

Sol(:,:) = [T1S_app(:), T1F_app(:), T2S_app(:), T2F_app(:), M0F_app(:), kFS_app(:), Delta_app(:)];
Bounds(:,:) = [T1S_LB(:), T1S_UB(:), T1F_LB(:), T1F_UB(:), T2S_LB(:), T2S_UB(:), T2F_LB(:), T2F_UB(:), M0F_LB(:), M0F_UB(:), kFS_LB(:), kFS_UB(:), Delta_LB(:), Delta_UB(:)];
Subsets(:,:) = [T1S_Sub(1:5), T1F_Sub(1:5), T2S_Sub(1:5), T2F_Sub(1:5), M0F_Sub(1:5), kFS_Sub(1:5), Delta_Sub(1:5)];

end

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