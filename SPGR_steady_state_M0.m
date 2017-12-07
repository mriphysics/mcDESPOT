% Using matrix exponential approach from Liu (2015).

function Mss_Sig = SPGR_steady_state_M0(FA, TR, varargin)

% Check for additional aprameters.
for ii = 1:length(varargin)
    
    if strcmpi(varargin{ii},'T1_S')
        T1_S = varargin{ii+1};
    end
    if strcmpi(varargin{ii},'T1_F')
        T1_F = varargin{ii+1};
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

% Calculate remaining parameters.
R1_S = 1/T1_S; R1_F = 1/T1_F;
%M0_S = 1 - M0_F; k_SF = (M0_F*k_FS)/M0_S; 

% Define Bloch-McConnell terms.
C = [(R1_F * M0_F) ; (R1_S * M0_S)];
A = [-R1_F-k_FS k_SF ; k_FS -R1_S-k_SF];

Mss = zeros(length(FA),2); Mss_Sig = zeros(length(FA),1);

AinvC = A\C;
em = expm(A * TR);

for jj = 1:length(FA)
    
    T = diag([cos(FA(jj)) cos(FA(jj))]);
    Mss(jj,:) = (eye(2) - (em * T))^-1 * (em - eye(2)) * AinvC;
    
    % Extract signal component.
    Mss_Sig(jj) = sin(FA(jj)) * (Mss(jj,1) + Mss(jj,2));
    
end

end