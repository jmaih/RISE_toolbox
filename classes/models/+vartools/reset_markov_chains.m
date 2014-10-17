function [markov_chains_]=reset_markov_chains(param_template,markov_chains_,endo_names)
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
% More About
% ------------
%
% Examples
% ---------
%
% See also: 


% create the parameters
%-----------------------
plist=vartools.create_baseline_parameters(param_template);

% add the transition probabilities
%---------------------------------
for im_=1:numel(markov_chains_)
    name=markov_chains_(im_).name;
    if strcmp(name,'const')
        error(['it is not allowed to define a markov chain with name "',const,'"'])
    end
    number_of_states=numel(markov_chains_(im_).states_expected_duration);
    for st=1:number_of_states
        for slead=1:number_of_states
            if st==slead
                continue
            end
            newparam=sprintf('%s_tp_%0.0f_%0.0f',name,st,slead);
            plist=[plist;newparam];
        end
    end
end

% push all the parameters into the markov chains
%------------------------------------------------
markov_chains_=vartools.format_markov_chains(markov_chains_,plist,endo_names);

end