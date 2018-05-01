function [lin_restr,nonlin_restr,markov_chains]=create_restrictions_and_markov_chains1(markov_chains)
% create_lin_restrictions_and_markov_chains1 -- creates restrictions and
% markov chains for the SVAR model in which coefficients are switching
% across all equations (synchronized case)
%
% ::
%
%
%   [lin_restr,nonlin_restr,tpl]=create_lin_restrictions_and_markov_chains1(tpl)
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

% We add a Markov chain to the template
%----------------------------------------
% N.B: The chain controls the coefficients of all equations but not the
% variance, which remains constant.
markov_chains(end+1)=struct('name','syncoef',...
    'states_expected_duration',[3+1i,3+1i],...
    'controlled_parameters',{{'c','a0','a1','a2'}});

% syntax is coef(eqtn,vname,lag,chain_name,state)
%------------------------------------------------
lin_restr=cell(0,1);
nonlin_restr=cell(0,1);

numberOfStates=numel(markov_chains(end).states_expected_duration);

for istate=1:numberOfStates
    
    mystate=int2str(istate);
    
    lin_restr=[lin_restr
        {
        % first equation or "FFR" equation:
        %----------------------------------
        ['a1(1,pi,syncoef,',mystate,')=0']
        ['a2(1,pi,syncoef,',mystate,')=0']
        ['a1(1,ygap,syncoef,',mystate,')=0']
        ['a2(1,ygap,syncoef,',mystate,')=0']
        ['a2(1,FFR,syncoef,',mystate,')=0']
        % second equation or "pi" equation
        %----------------------------------
        ['a0(2,FFR,syncoef,',mystate,')=0']
        ['a1(2,FFR,syncoef,',mystate,')=0']
        ['a2(2,FFR,syncoef,',mystate,')=0']
        ['a1(2,ygap,syncoef,',mystate,')=0']
        ['a2(2,ygap,syncoef,',mystate,')=0']
        % third equation or "ygap" equation
        %-----------------------------------
        ['a1(3,FFR,syncoef,',mystate,')=0']
        ['a2(3,FFR,syncoef,',mystate,')=0']
        ['a1(3,pi,syncoef,',mystate,')=0']
        ['a2(3,pi,syncoef,',mystate,')=0']
        ['a0(3,pi,syncoef,',mystate,')+a0(3,FFR,syncoef,',mystate,')=0']
        }
        ]; %#ok<AGROW>
    nonlin_restr=[nonlin_restr
        {
        ['a0(3,FFR,syncoef,',mystate,')>=0']
        ['a1(1,FFR,syncoef,',mystate,')>=0']
        ['a1(1,FFR,syncoef,',mystate,')<=1']
        }
        ]; %#ok<AGROW>
end

end