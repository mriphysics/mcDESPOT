% Using matrix exponential approach from Liu (2015).

function Mss_Sig = CS_SPGR_steady_state_MT_M0(FA, TR, varargin)

% Calculate bound pool saturation factor.
G = 1.4e-5; % [Gloor]
Gamma = 2 * pi * 42.57747892e6; % [rad/s/T]
B1 = 13e-6; % [T]
T_RF = deg2rad(50)/(Gamma * B1);
W = (pi/T_RF) * (Gamma * B1)^2 * T_RF * G;

% Check for additional arguments.
for ii = 1:length(varargin)
    
    if strcmpi(varargin{ii},'T1_S')
        T1_S = varargin{ii+1};
    end
    if strcmpi(varargin{ii},'T1_F')
        T1_F = varargin{ii+1};
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
k_SF = (M0_F * k_FS)/M0_S ; k_BS = (M0_S * k_SB)/M0_B ; k_BF = (M0_F * k_FB)/M0_B ;
R1_S = 1/T1_S; R1_F = 1/T1_F; R1_B = 1/T1_B;

% Define Bloch-McConnell terms.
C = [(R1_F * M0_F) ; (R1_S * M0_S) ; (R1_B * M0_B)];
A = [-R1_F-k_FS-k_FB k_SF k_BF; k_FS -R1_S-k_SF-k_SB k_BS ; k_FB k_SB -R1_B-k_BF-k_BS];

Mss = zeros(length(FA),3);
Mss_Sig = zeros(length(FA),1);

for jj = 1:length(FA)

    AinvC = A\C;
    em = expm(A * TR);
    T = [cos(FA(jj)) 0 0 ; 0 cos(FA(jj)) 0 ; 0 0 exp(-W * T_RF)]; 
    Mss(jj,:) = (eye(3) - (em * T))^-1 * (em - eye(3)) * AinvC; 
    % Extract signal component.
    Mss_Sig(jj) = sin(FA(jj)) * (Mss(jj,1) + Mss(jj,2));
    
end

% dMzB = -k_FB*Mss(jj,1) + k_SB*Mss(jj,2) - (R1_B + k_BF + k_BS + W)*Mss(jj,3) + M0_B*R1_B;

end
