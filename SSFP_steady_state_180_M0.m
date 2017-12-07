% Using matrix exponential approach from Liu (2015).

function Mss_Sig = SSFP_steady_state_180_M0(FA, TR, varargin)

% Check for these additional arguments specified in separate script.
for ii = 1:length(varargin)
    
    if strcmpi(varargin{ii},'T1_S')
        T1_S = varargin{ii+1};
    end
    if strcmpi(varargin{ii},'T2_S')
        T2_S = varargin{ii+1};
    end
    if strcmpi(varargin{ii},'T1_F')
        T1_F = varargin{ii+1};
    end
    if strcmpi(varargin{ii},'T2_F')
        T2_F = varargin{ii+1};
    end
    if strcmpi(varargin{ii},'M0_F')
        M0_F = varargin{ii+1};
    end
    if strcmpi(varargin{ii},'M0_S')
      M0_S = varargin{ii+1};
    end
    if strcmpi(varargin{ii},'k_FS')
        k_FS = varargin{ii+1};
    end
    if strcmpi(varargin{ii},'k_SF')
      k_SF = varargin{ii+1};
    end
end

% Inclusion of relaxation. Define relaxation rates for simplification.
R1_S = 1/T1_S; R2_S = 1/T2_S; R1_F = 1/T1_F; R2_F = 1/T2_F;
%M0_S = 1 - M0_F; k_SF = (M0_F*k_FS)/M0_S;
C = [0 ; 0 ; (M0_F * R1_F) ; 0 ; 0 ; (M0_S * R1_S)];

A = [-R2_F-k_FS 0 0 k_SF 0 0 ; 0 -R2_F-k_FS 0 0 k_SF 0 ;
    0 0 -R1_F-k_FS 0 0 k_SF ; k_FS 0 0 -R2_S-k_SF 0 0 ;
    0 k_FS 0 0 -R2_S-k_SF 0 ; 0 0 k_FS 0 0 -R1_S-k_SF];

% Accounts for phase-cycling.
PC_M = [-1 0 0 0 0 0 ; 0 -1 0 0 0 0 ; 0 0 1 0 0 0 ; 0 0 0 -1 0 0 ; 0 0 0 0 -1 0 ; 0 0 0 0 0 1];

Mss = zeros(length(FA),6);
Mss_Sig = zeros(length(FA),1);

AinvC = A\C;
em = expm(A * TR);

for jj = 1:length(FA)
    
    T = [1 0 0 0 0 0 ; 0 cos(FA(jj)) sin(FA(jj)) 0 0 0 ; 0 -sin(FA(jj)) cos(FA(jj)) 0 0 0 ;
        0 0 0 1 0 0 ; 0 0 0 0 cos(FA(jj)) sin(FA(jj)) ; 0 0 0 0 -sin(FA(jj)) cos(FA(jj))];
    
    Mss(jj,:) = (PC_M - (em * T))^-1 * (em - eye(6)) * AinvC;
    
    Mss_Sig(jj) = abs((Mss(jj,1) + (1i * Mss(jj,2))) + (Mss(jj,4) + (1i * Mss(jj,5))));
    
end

end