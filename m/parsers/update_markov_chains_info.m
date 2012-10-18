function [MC,istp]=update_markov_chains_info(MC,parameter,flag)

% creates a markov chain info and updates it as new information comes in.
if nargin<3
    flag=false; % endogenous or exogenous switching
    if nargin<2
        if nargin<1
            MC=[];
        end
    end
end
% order is: 1)chain name, 2) # regimes, 3) status (endo,exo), 4) script
if isempty(MC)
    MC=transpose({'const',1,0,'Qi=1;'}); 
    % constants have at most one state cell(2,0);
    % determine the number of markov chains using the parameter list... and not
    % ParameterArray. Then check that each parameter in ParameterArray is
    % assigned a regime
end
istp=false;
if nargin>1
    [istp,junk,chain_name,max_state]=is_transition_probability(parameter);
    if istp
        loc=find(strcmp(chain_name,MC(1,:)));
        if isempty(loc)
            MC=[MC,transpose({chain_name,max_state,flag,''})];
        else
            MC{2,loc}=max(max_state,MC{2,loc});
            if flag~=MC{3,loc}
                disp([mfilename,':: it is forbidden to list endogenous chains as parameters...'])
                disp('Use auxiliary parameters if one of the probabilities is constant')
                error([mfilename,':: Markov chain ''',chain_name,''' cannot be both endogenous and exogenous'])
            end
        end
    end
end

