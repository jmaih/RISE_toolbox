function obj=add_markov_chains_and_parameters(obj,markov_chains_)
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
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

if ~isfield(markov_chains_,'param_list')
    return
end
plist={};
plist_tex={};
gov_chain=[];
is_switching=[];
for ichain=1:numel(markov_chains_)
    thisList=markov_chains_(ichain).param_list(:)';
    thisList_tex=markov_chains_(ichain).param_list_tex(:)';
    nlist=numel(thisList);
    plist=[plist,thisList]; %#ok<*AGROW>
    plist_tex=[plist_tex,thisList_tex];
    gov_chain=[gov_chain,ichain*ones(1,nlist)];
    is_switching=[is_switching,markov_chains_(ichain).is_switching*ones(1,nlist)];
end
[obj.parameters.name,tags]=sort(plist);
obj.parameters.governing_chain=gov_chain(tags);
obj.parameters.is_switching=logical(is_switching(tags));
obj.parameters.tex_name=plist_tex(tags);
% re-tag the transition probabilities
%------------------------------------
obj.parameters.number=numel(obj.parameters.name);
for iparam=1:obj.parameters.number
    obj.parameters.is_trans_prob(1,iparam)=...
        parser.is_transition_probability(obj.parameters.name{iparam});
end
[~,obj.routines.transition_matrix,obj.markov_chains]=...
    parser.transition_probabilities({'param'},obj.parameters.name,markov_chains_,'','',[]);
end

