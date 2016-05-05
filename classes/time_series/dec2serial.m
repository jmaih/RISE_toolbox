function s=dec2serial(varargin)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

narginchk(1,3)

if nargin==1
    
    year=[varargin{1}.year];
    
    period=[varargin{1}.period];
    
    freq=[varargin{1}.freq];
    
else
    
    year=varargin{1};
    
    period=varargin{2};
    
    freq=varargin{3};
    
end

stamp=time_frequency_stamp();

dn=@(x,per,freq)x.*freq+per-1+stamp(freq);

s=dn(year,period,freq);

end