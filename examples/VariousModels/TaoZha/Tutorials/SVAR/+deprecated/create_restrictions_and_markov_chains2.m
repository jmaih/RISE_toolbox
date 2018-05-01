function [lin_restr,nonlin_restr,markov_chains]=create_restrictions_and_markov_chains2(markov_chains)
% create_restrictions_and_markov_chains2 -- creates restrictions and
% markov chains for the SVAR model in which Coefficients and variances have
% different chains, different regimes, and different durations
%
% ::
%
%
%   [lin_restr,nonlin_restr,tpl]=create_restrictions_and_markov_chains2(tpl)
%
% Args:
%
%    - **markov_chains** [empty|struct]: structure of previously defined
%    markov chains
%
% Returns:
%    :
%
%    - **lin_restr** [cell]: cell array of restrictions (see below).
%
%    - **nonlin_restr** [cell]: cell array of inequality restrictions
%
%    - **markov_chains** [struct]: modified markov chains
%
% Note:
%
%    - The syntax to construct a restriction
%      --> ai(eqtn)
%      --> ai(eqtn,vbl)
%      --> ai(eqtn,vbl,chain_name,state)
%      --> a(eqtn)
%      --> a(eqtn,vbl)
%      --> a(eqtn,vbl,chain_name,state)
%      - **eqtn** [integer]: integer
%      - **vbl** [integer|char]: integer or variable name
%      - **i** [integer]: lag
%      - **chain_name** [char]: name of the markov chain
%      - **state** [integer]: state number
%
%    - The lag coefficients are labelled a0, a1, a2,...,ak, for a model with k
%    lags. Obviously, a0 denotes the contemporaneous coefficients.
%
%    - The constant terms labelled c_1_1, c_2_2,...,c_n_n, for a model with n
%    endogenous variables.
%
%    - The standard deviations labelled s_1_1, s_2_2,...,s_n_n, for a
%    model with n endogenous variables.
%
% Example:
%
%    See also:

if nargin==0||isempty(markov_chains)
    
    markov_chains=struct('name',{},...
    'states_expected_duration',{},...
    'controlled_parameters',{});
    
end

% We borrow both the restrictions and the markov chains from the model in
% which all coefficients across all equations switch in lockstep.
%--------------------------------------------------------------------------
[lin_restr,nonlin_restr,markov_chains]=create_restrictions_and_markov_chains1(markov_chains);

% Then we add another chain controling all variances across all equations
%-------------------------------------------------------------------------
% We just make sure we do not change the restrictions
last=numel(markov_chains);

markov_chains(last+1)=struct('name','syncvol',...
    'states_expected_duration',[2+1i,2+1i,2+1i],...
    'controlled_parameters',{{'s'}});

end