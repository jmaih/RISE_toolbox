function Linf=find_long_run(B,L0,nlags,endo_nbr)
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

beta=0;
for ilag = 1:nlags
    beta_i=B(:,(ilag-1)*endo_nbr+1:ilag*endo_nbr);
    beta = beta + beta_i;
end
Linf = (eye(endo_nbr)-beta)\L0;
end
