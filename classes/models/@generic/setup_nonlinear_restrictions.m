function obj=setup_nonlinear_restrictions(obj)
% setup_nonlinear_restrictions - sets nonlinear restrictions
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
%    - uses estim_nonlinear_restrictions, which should be a cell array. Each
%    item of the array is a string of the form
%      - 'f(p1,p2,...,pn)>=h(p1,p2,...,pn)'
%      - 'f(p1,p2,...,pn)>h(p1,p2,...,pn)'
%      - 'f(p1,p2,...,pn)<=h(p1,p2,...,pn)'
%      - 'f(p1,p2,...,pn)<h(p1,p2,...,pn)'
%      - 'pj=h(p1,p2,...,pn)'
%
%    - In the statements above,
%      - eqtn [digits|variable name]
%      - vbl [digits|variable name]
%      - lag [digits]
%      - chain [char]
%      - state [digits]
%
% Example:
%
%    See also:

RestrictionsBlock=obj.options.estim_nonlinear_restrictions;

if isstruct(RestrictionsBlock)
    
    RestrictionsBlock=RestrictionsBlock.original;
    
end

if ~isempty(RestrictionsBlock)
    
    RestrictionsBlock=cellfun(@(x)x(~isspace(x)),...
        RestrictionsBlock,'uniformOutput',false);
    
end

param_names=obj.parameters.name;

governing_chain=obj.parameters.governing_chain;

chain_names=obj.markov_chains.small_markov_chain_info.chain_names;

regimes=cell2mat(obj.markov_chains.small_markov_chain_info.regimes(2:end,2:end));

[nonlinear_restrictions,is_linear_restriction,derived_parameters]=...
    generic_tools.nonlinear_restrictions_engine(...
    param_names,regimes,chain_names,governing_chain,...
    RestrictionsBlock);

obj.number_of_restrictions=struct('auxiliary',sum(is_linear_restriction),...
    'linear',nan,...
    'nonlinear',sum(~is_linear_restriction));

obj=add_to_routines(obj,'derived_parameters',derived_parameters,...
    'nonlinear_restrictions',nonlinear_restrictions);

end