% Using matrix exponential approach from Liu (2015).

function Mss_Sig = TRF_SSFP_steady_state_MT_M0(FA, TR, varargin)

% Calculate bound pool saturation factor.
G = 1.4e-5; % [Gloor]
Gamma = 2 * pi * 42.57747892e6; % [rad/s/T]
B1 = 13e-6; % [T]

% Check for additional arguments.
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
    if strcmpi(varargin{ii},'T1_B')
        T1_B = varargin{ii+1};
    end
    if strcmpi(varargin{ii},'M0_B')
        M0_B = varargin{ii+1};
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
    if strcmpi(varargin{ii},'k_SB')
        k_SB = varargin{ii+1};
    end
    if strcmpi(varargin{ii},'k_FB')
        k_FB = varargin{ii+1};
    end
    
end

% Calculate remaining parameters.
k_SF = (M0_F * k_FS)/M0_S; k_BS = (M0_S * k_SB)/M0_B ; k_BF = (M0_F * k_FB)/M0_B ;
R1_S = 1/T1_S; R2_S = 1/T2_S; R1_F = 1/T1_F; R2_F = 1/T2_F; R1_B = 1/T1_B;

% Define Bloch-McConnell terms.
C = [0 ; 0 ; (M0_F * R1_F) ; 0 ; 0 ; (M0_S * R1_S) ; (M0_B * R1_B)];

A = [-R2_F-k_FS 0 0 k_SF 0 0 0 ; 0 -R2_F-k_FS 0 0 k_SF 0 0 ;
    0 0 -R1_F-k_FS-k_FB 0 0 k_SF k_BF; k_FS 0 0 -R2_S-k_SF 0 0 0 ;
    0 k_FS 0 0 -R2_S-k_SF 0 0 ; 0 0 k_FS 0 0 -R1_S-k_SF-k_SB k_BS ; 0 0 k_FB 0 0 k_SB -R1_B-k_BF-k_BS];

Mss = zeros(length(FA),7);
Mss_Sig = zeros(length(FA),1);

AinvC = A\C;
em = expm(A * TR);

for jj = 1:length(FA)
    
    T_RF = FA(jj)/(Gamma * B1);
    W = (pi/T_RF) * (Gamma * B1)^2 * T_RF * G;
    
    T = [1 0 0 0 0 0 0 ; 0 cos(FA(jj)) sin(FA(jj)) 0 0 0 0 ; 0 -sin(FA(jj)) cos(FA(jj)) 0 0 0 0 ;
        0 0 0 1 0 0 0 ; 0 0 0 0 cos(FA(jj)) sin(FA(jj)) 0 ; 0 0 0 0 -sin(FA(jj)) cos(FA(jj)) 0 ; 0 0 0 0 0 0 exp(-W * T_RF)];
    
    Mss(jj,:) = ((eye(7) - (em * T))^-1) * (em - eye(7)) * AinvC;
    
    Mss_Sig(jj) = abs((Mss(jj,1) + (1i * Mss(jj,2))) + (Mss(jj,4) + (1i * Mss(jj,5))));
    
end

end