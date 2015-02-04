function flag=exitflag(flag,x1,fx1,TolFun)
% exitflag -- recast the exit flag from fsolve and lsqnonlin
%
% Syntax
% -------
% ::
%
%   flag=exitflag(flag,x1,fx1,TolFun)
%
% Inputs
% -------
%
% - **flag** [scalar]: exit flag from fsolve or lsqnonlin
%
% - **x1** [scalar|vector|matrix]: optimum of the function to set to zero
%   (f(x)=0)
%
% - **fx1** [scalar]: value or norm of f(x1)
%
% - **TolFun** [{sqrt(eps)}|scalar]: scalar or vector of model objects
%
% Outputs
% --------
%
% - **flag** [scalar]: recomputed exitflag
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

if nargin<4
    TolFun=sqrt(eps);
end

if ismember(flag,[1:4,-3]) && utils.error.valid(x1)
    flag=inf;
    if fx1 <=TolFun
        flag=1;
        %         x1=reshape(x1,size(x0));
    end
end
