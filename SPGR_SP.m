function Mss_Sig = SPGR_SP(FA, TR, varargin)

for ii = 1:length(varargin)
    
    if strcmpi(varargin{ii},'T1') 
        T1 = varargin{ii+1};  
    end

    if strcmpi(varargin{ii},'M0') 
        M0 = varargin{ii+1};  
    end
    
end

R1 = 1/T1;

C = (R1 * M0);
A = -R1;

Mss = zeros(length(FA),1); Mss_Sig = zeros(length(FA),1);

AinvC = A\C;
em = expm(A * TR);

for jj = 1:length(FA)
    
    T = cos(FA(jj));
    
    Mss(jj) = (eye(1) - (em * T))^-1 * (em - eye(1)) * AinvC;
    
    Mss_Sig(jj) = sin(FA(jj)) * Mss(jj);
    
end

end