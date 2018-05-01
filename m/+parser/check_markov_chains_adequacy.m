function check_markov_chains_adequacy(dictionary)
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


% check that the markov chain is parameterized in all of its elements
%---------------------------------------------------------------------

markov_chain_names={dictionary.markov_chains.name};
tp_endo=dictionary.time_varying_probabilities;
param_names={dictionary.parameters.name};
tp_exo=param_names([dictionary.parameters.is_trans_prob]);
union_endo_exo=union(tp_endo,tp_exo);

% intersection of the markov chains must be empty
%------------------------------------------------

missing_trans_probs={};
for ichain=1:numel(markov_chain_names)
    if strcmp(markov_chain_names{ichain},'const')
        continue
    end
    nstates=dictionary.markov_chains(ichain).number_of_states;
    for istate=1:nstates
        for jstate=1:nstates
            if jstate~=istate
                this_prob=[markov_chain_names{ichain},...
                    '_tp_',sprintf('%0.0f',istate),'_',sprintf('%0.0f',jstate)];
                if ~ismember(this_prob,union_endo_exo)
                    missing_trans_probs=[missing_trans_probs,this_prob]; %#ok<AGROW>
                end
            end
        end
    end
end
if ~isempty(missing_trans_probs)
    disp(missing_trans_probs(:)')
    error('The transition probabilities above are missing')
end
