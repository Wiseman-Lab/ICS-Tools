%
function out = ICSFit(params,xiLags,etaLags,k_on,k_off,varargin)

errBool = 0;
includeZeroLag = 0;

for i = 1:length(varargin)
    if any(strcmpi(varargin{i},{'error','err','residual','res'}))
        errBool = 1;
        corrData = varargin{i+1};
    elseif any(strcmpi(varargin{i},{'includeZeroLag','zeroLag'}))
        includeZeroLag = 1;
    end
end

if errBool && ~includeZeroLag
    zeroLagInd = find(xiLags==0 & etaLags==0);
    xiLags(zeroLagInd) = [];
    etaLags(zeroLagInd) = [];
   
    corrData(zeroLagInd) = [];    
end

s = struct('density',params(1),'offset',params(2),'w1',params(3),'w2',params(3),'theta',0);

fields = fieldnames(s);

for n = 4:length(fields)
    if n <= length(params)   
        fieldname = fields{n};
        s.(fieldname) = params(n);
    end
end

K = k_on + k_off;
    
F = K/k_on*1/(s.density*pi*min(s.w1,s.w2)^2).*exp((-(xiLags.^2+etaLags.^2).*(s.w1^2+s.w2^2)+(s.w1^2-s.w2.^2).*((xiLags.^2-etaLags.^2)*cos(2*s.theta)-etaLags.*xiLags*sin(2*s.theta)))/(2*s.w1^2*s.w2^2))+s.offset;

% K/k_on*1/(s.density*pi*min(s.w1,s.w2)^2)

if errBool
    err = norm(corrData - F,2);
end

if ~errBool
    out = F;
else
    out = err;
end
