function [restrictions,markov_chains,switch_prior]=create_restrictions_and_markov_chains3(markov_chains,switch_prior)
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
        'number_of_states',{},...
        'controlled_parameters',{},...
        'endogenous_probabilities',{},...
        'probability_parameters',{});
    
    switch_prior=struct();
    
end

% We add a Markov chain to the template
%----------------------------------------
% N.B: The chain controls the coefficients of the policy (RRF) equation
% only. It does not control the variance
last=numel(markov_chains);

markov_chains(last+1)=struct('name','mpcoef',...
    'number_of_states',2,...
    'controlled_parameters',{{'c(1)','a0(1)','a1(1)','a2(1)'}},...
    'endogenous_probabilities',[],...
    'probability_parameters',[]);

switch_prior.mpcoef_tp_1_2={0.5,0.1,0.3,'beta'};
switch_prior.mpcoef_tp_2_1={0.5,0.1,0.3,'beta'};

% syntax is alag(eqtn,vname,chain_name,state)
%------------------------------------------------
lin_restr=cell(0,1);

nonlin_restr={'a0(3,FFR)>=0'};

% first equation or "FFR" equation: 
%----------------------------------
for istate=1:markov_chains(end).number_of_states
    
    mystate=int2str(istate);
    
    lin_restr=[lin_restr
        {
        ['a1(1,pi,mpcoef,',mystate,')=0']
        ['a2(1,pi,mpcoef,',mystate,')=0']
        ['a1(1,ygap,mpcoef,',mystate,')=0']
        ['a2(1,ygap,mpcoef,',mystate,')=0']
        ['a2(1,FFR,mpcoef,',mystate,')=0']
        }
        ]; %#ok<AGROW>
    
    nonlin_restr=[nonlin_restr
        {
        ['a1(1,FFR,mpcoef,',mystate,')>=0']
        ['a1(1,FFR,mpcoef,',mystate,')<=1']
        }
        ]; %#ok<AGROW>
end

lin_restr=[lin_restr
    {
    % second equation or "pi" equation
    %----------------------------------
    'a0(2,FFR)=0'
    'a1(2,FFR)=0'
    'a2(2,FFR)=0'
    'a1(2,ygap)=0'
    'a2(2,ygap)=0'
    % third equation or "ygap" equation
    %-----------------------------------
    'a1(3,FFR)=0'
    'a2(3,FFR)=0'
    'a1(3,pi)=0'
    'a2(3,pi)=0'
    'a0(3,pi)+a0(3,FFR)=0'
    }
    ];

restrictions=[lin_restr;nonlin_restr];

end