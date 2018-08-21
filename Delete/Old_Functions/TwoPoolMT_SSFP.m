% Using matrix exponential approach from Liu (2015).

function Mss_Sig = TwoPoolMT_SSFP(FA, TR, varargin)

% Calculate bound pool saturation factor.
G = 1.4e-5; % [Gloor]
Gamma = 2 * pi * 42.57747892e6; % [rad/s/T]
B1 = 13e-6; % [T]

T_RF = deg2rad(50)/(Gamma * B1);
W = (pi/T_RF) * (Gamma * B1)^2 * T_RF * G;

% T_RF = max(FA)/(Gamma * B1);

% Check for additional arguments.
for ii = 1:length(varargin)
    
    if strcmpi(varargin{ii},'T1_W')
        T1_W = varargin{ii+1};
    end
    if strcmpi(varargin{ii},'T2_W')
        T2_W = varargin{ii+1};
    end
    if strcmpi(varargin{ii},'T1_B')
        T1_B = varargin{ii+1};
    end
    if strcmpi(varargin{ii},'M0_B')
        M0_B = varargin{ii+1};
    end
    if strcmpi(varargin{ii},'k_WB')
        k_WB = varargin{ii+1};
    end
    
end

% Calculate remaining parameters.
M0_W = 1 - M0_B; k_BW = (M0_W * k_WB)/M0_B ;
R1_W = 1/T1_W; R2_W = 1/T2_W; R1_B = 1/T1_B;

% Define Bloch-McConnell terms.
C = [0 ; 0 ; (R1_W*M0_W) ; (R1_B*M0_B)];

A = [-R2_W 0 0 0 ; 0 -R2_W 0 0 ; 0 0 -(R1_W + k_WB) k_BW ; 0 0 k_WB -(R1_B + k_BW)];

Mss = zeros(length(FA),4);
Mss_Sig = zeros(length(FA),1);

AinvC = A\C;
em = expm(A * TR);

for jj = 1:length(FA)
    
%     T_RF = FA(jj)/(Gamma * B1);
%     W = (pi/T_RF) * (Gamma * B1)^2 * T_RF * G;

%     SF = FA(jj)/max(FA);
%     W = (pi/T_RF) * (Gamma * B1 * SF)^2 * T_RF * G;
    
    T = [1 0 0 0 ; 0 cos(FA(jj)) sin(FA(jj)) 0 ; 0 -sin(FA(jj)) cos(FA(jj)) 0 ; 0 0 0 exp(-W * T_RF)];
    
    Mss(jj,:) = ((eye(4) - (em * T))^-1) * (em - eye(4)) * AinvC;
    
    Mss_Sig(jj) = abs((Mss(jj,1) + (1i * Mss(jj,2))));
    
end

end