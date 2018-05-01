function shocks=set_exogenous_data(exo_nbr,is_observed,simul_with_shocks,z,horizon)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

if nargin<5
    horizon=1;
end
shocks=zeros(exo_nbr,horizon);
if simul_with_shocks
    shocks(~is_observed,:)=randn(sum(~is_observed),horizon);
end
% deterministic variables
%------------------------
ncols=size(z,2);
if ncols
    cutz=min(ncols,horizon);
    shocks(is_observed,1:cutz)=z(:,1:cutz);
end
