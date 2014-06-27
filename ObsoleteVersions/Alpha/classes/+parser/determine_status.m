function [status,loc_]=determine_status(dictionary,tokk)
loc_=[];
possibilities={'orig_endogenous','y'
    'exogenous','x'
    'parameters','param' % same as definitions
    'definitions','def'
    'time_varying_probabilities','tvp'
    'known_words','f' % same as functions
    'symbols',tokk
    'add_operators','+'
    'mult_operators','*'
    'relational_operators','>'
    'chain_names','cn'
    };
iter=0;
status='unknown';
while isempty(loc_) && iter<size(possibilities,1)
    iter=iter+1;
    if iter<4
        loc_=find(strcmp(tokk,{dictionary.(possibilities{iter,1}).name}),1);
    else
        loc_=find(strcmp(tokk,dictionary.(possibilities{iter,1})),1);
    end
    if ~isempty(loc_)
        status=possibilities{iter,2};
        break
    end
end
if isempty(loc_)
    if exist([tokk,'.m'],'file')
        status='f'; % same as known words
    else
        try %#ok<TRYNC>
            isnumeric(eval(tokk));
            status='n';
        end
    end
    % maybe it is a steady state definition
    if strncmp(tokk,'xx_ssmdef_',10) && all(isstrprop(tokk(11:end),'digit'))
        status='f'; % same as known words and functions
    end
end
% maybe this function should be called with a block name...
% %         if strcmp(status,'unknown') && strcmp(current_block_name,'steady_state_model')
% %             dictionary.known_words=[dictionary.known_words,{tokk}];
% %             status='f';
% %         end
end
