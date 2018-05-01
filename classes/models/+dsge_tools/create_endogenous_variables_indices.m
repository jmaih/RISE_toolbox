function [the_leads,the_lags,nind,endo_nbr]=create_endogenous_variables_indices(lead_lag_incidence)

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

endo_nbr=size(lead_lag_incidence,1);

the_leads=find(lead_lag_incidence(:,1)>0);

the_current=(1:endo_nbr)';

the_lags=find(lead_lag_incidence(:,3)>0);

indices=[the_leads;the_current;the_lags];

nind=numel(indices);

end