function [lin_restr,nonlin_restr,tpl]=create_restrictions_and_markov_chains3(tpl)
% create_restrictions_and_markov_chains3 -- creates restrictions and
% markov chains for the SVAR model in which only the parameters for the
% monetary policy equation are changing.
%
% ::
%
%
%   [lin_restr,tpl]=create_restrictions_and_markov_chains3(tpl)
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

% We add a Markov chain to the template
%----------------------------------------
% N.B: The chain controls the coefficients of the policy (RRF) equation
% only. It does not control the variance
last=numel(tpl.markov_chains);
tpl.markov_chains(last+1)=struct('name','mpcoef',...
    'states_expected_duration',[3+1i,3+1i],...
    'controlled_parameters',{{'c(1)','a0(1)','a1(1)','a2(1)'}});

% syntax is coef(eqtn,vname,lag,chain_name,state)
%------------------------------------------------
lin_restr=cell(0,2);

numberOfStates=numel(tpl.markov_chains(end).states_expected_duration);

% first equation or "FFR" equation: 
%----------------------------------
for istate=1:numberOfStates
    lin_restr=[lin_restr
        {
        coef('FFR','pi',1,'mpcoef',istate),0
        coef('FFR','pi',2,'mpcoef',istate),0
        coef('FFR','ygap',1,'mpcoef',istate),0
        coef('FFR','ygap',2,'mpcoef',istate),0
        coef('FFR','FFR',2,'mpcoef',istate),0
        }
        ]; %#ok<AGROW>
end
lin_restr=[lin_restr
    {
    % second equation or "pi" equation
    %----------------------------------
    coef('pi','FFR',0),0
    coef('pi','FFR',1),0
    coef('pi','FFR',2),0
    coef('pi','ygap',1),0
    coef('pi','ygap',2),0
    % third equation or "ygap" equation
    %-----------------------------------
    coef('ygap','FFR',1),0
    coef('ygap','FFR',2),0
    coef('ygap','pi',1),0
    coef('ygap','pi',2),0
    coef('ygap','pi',0)+coef('ygap','FFR',0),0
    }
    ];

% borrow nonlinear restrictions from the constant-parameter model
% since they only concern the third equation in which parameters
% do not switch
%----------------------------------------------------------------
% do not overwrite the linear restrictions or the template
[~,nonlin_restr]=create_restrictions_and_markov_chains0(tpl);
end