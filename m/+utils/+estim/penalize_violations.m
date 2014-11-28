function p=penalize_violations(viol,c)
% penalize_violations - penalizes the violation of general restrictions
%
% Syntax
% -------
% ::
%
%   p=penalize_violations(viol)
%   p=penalize_violations(viol,c)
%
% Inputs
% -------
%
% - **viol** [vector] : value of constraints expressed in the form g(x)<=0
%
% - **c** [scalar|{10}]: penalty factor
%
% Outputs
% --------
%
% - **p** [scalar] : penalty
%
% More About
% ------------
%
% - restrictions are expected to be of the form g(x)<=0
%
% Examples
% ---------
%
% See also:

if nargin<2
    c=[];
end
if isempty(c)
    c=10;
end

if c<=0
    error('penalty factor must be strictly positive')
end

p=0;
if ~isempty(viol)
    viol=max(0,viol(:));
    viol=viol(viol>0);
    if ~isempty(viol)
        p=c*sum((viol/min(viol)).^2);
    end
end

end