function [lin_restr,nonlin_restr,tpl]=create_restrictions_and_markov_chains6(tpl)
% create_restrictions_and_markov_chains6 -- creates linear restrictions and
% markov chains for the SVAR model in which only variances in all 3
% equations are switching.
%
% ::
%
%
%   [lin_restr,nonlin_restr,tpl]=create_restrictions_and_markov_chains6(tpl)
%
% Args:
%
%    - **tpl** [struct]: template created for SVAR objects
%
% Returns:
%    :
%
%    - **lin_restr** [cell]: two column-cell (see below). The first column
%    contains COEF objects or linear combinations of COEF objects, which are
%    themselves COEF objects.
%
%    - **nonlin_restr** [cell]: one column-cell containing inequality or
%    nonlinear restrictions
%
%    - **tpl** [struct]: modified template
%
% Note:
%
%    - The syntax to construct an advanced COEF object is
%    a=coef(eqtn,vname,lag,chain_name,state)
%      - **eqtn** [integer|char]: integer or variable name
%      - **vname** [integer|char]: integer or variable name
%      - **lag** [integer]: integer or variable name
%      - **chain_name** [char]: name of the markov chain
%      - **state** [integer]: state number
%
%    - RISE sorts the endogenous variables alphabetically and use this order
%    to tag each equation in SVAR and RFVAR models.
%
%    - The lag coefficients are labelled a0, a1, a2,...,ak, for a model with k
%    lags. Obviously, a0 denotes the contemporaneous coefficients.
%
%    - The constant terms labelled c_1_1, c_2_2,...,c_n_n, for a model with n
%    endogenous variables.
%
%    - The standard deviations labelled sig_1_1, sig_2_2,...,sig_n_n, for a
%    model with n endogenous variables.
%
% Example:
%    coef('pi','ygap',0,'policy',1)
%
%    coef(2,3,0,'policy',1)
%
%    coef(2,'ygap',0,'policy',1)
%
%    coef('pi',3,0,'policy',1)
%
%    See also:

% We borrow both the restrictions and the markov chains from the model in
% which coefficients do not switch.
%--------------------------------------------------------------------------
[lin_restr,nonlin_restr,tpl]=create_restrictions_and_markov_chains0(tpl);

% Then we add a chain controling all variances across all equations
%--------------------------------------------------------------------
% We just make sure we do not change the restrictions
last=numel(tpl.markov_chains);
tpl.markov_chains(last+1)=struct('name','syncvol',...
    'states_expected_duration',[2+1i,2+1i,2+1i],...
    'controlled_parameters',{{'sig'}});

end