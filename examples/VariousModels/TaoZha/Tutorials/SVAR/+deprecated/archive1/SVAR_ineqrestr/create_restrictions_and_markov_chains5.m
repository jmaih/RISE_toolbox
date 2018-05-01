function [lin_restr,nonlin_restr,tpl]=create_restrictions_and_markov_chains5(tpl)
% create_restrictions_and_markov_chains5 -- creates restrictions and
% markov chains for the SVAR model in which both coefficients and variance
% for the monetary policy equation are changing with two independent Markov
% processes
%
% ::
%
%
%   [lin_restr,nonlin_restr,tpl]=create_restrictions_and_markov_chains5(tpl)
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

% The parameter restrictions are identical to those in the model
% with regime switching in the policy coefficients only. Hence the mp_coef
% markov chain will also be common to those two models.
[lin_restr,nonlin_restr,tpl]=create_restrictions_and_markov_chains3(tpl);

% We add the volatility Markov chain from the model in which only the
% volatility of the monetary policy equation changes. N.B: In the process,
% we want to make sure we do not over-write the restrictions above!
%--------------------------------------------------------------------------
[~,~,tpl]=create_restrictions_and_markov_chains4(tpl);

end