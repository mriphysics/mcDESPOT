function Mss_Sig = SSFP_SP_180(FA, TR, varargin)

for ii = 1:length(varargin)
    
    if strcmpi(varargin{ii},'T1')
        T1 = varargin{ii+1};
    end
    
    if strcmpi(varargin{ii},'T2')
        T2 = varargin{ii+1};
    end
    
    if strcmpi(varargin{ii},'M0')
        M0 = varargin{ii+1};
    end
    
end

R1 = 1/T1; R2 = 1/T2; 

C = [0 ; 0 ; (M0 * R1)];
A = [-R2 0 0 ; 0 -R2 0 ; 0 0 -R1];

PC_M = [-1 0 0 ; 0 -1 0 ; 0 0 1];

Mss = zeros(length(FA),3);
Mss_Sig = zeros(length(FA),1);

AinvC = A\C;
em = expm(A * TR);

for jj = 1:length(FA)
    
    T = [1 0 0 ; 0 cos(FA(jj)) sin(FA(jj)) ; 0 -sin(FA(jj)) cos(FA(jj))];
    
    Mss(jj,:) = (PC_M - (em * T))^-1 * (em - eye(3)) * AinvC;
    
    Mss_Sig(jj) = abs((Mss(jj,1) + (1i * Mss(jj,2))));
    
end

end