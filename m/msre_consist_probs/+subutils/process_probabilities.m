function inner_probs=process_probabilities(obj,inner_probs)

if all(cellfun(@isstruct,inner_probs(:,2),'uniformOutput',true))
    
    return
    
end

regimes=cell2mat(obj.markov_chains.regimes(2:end,2:end));

ordered_names=obj.endogenous.name(obj.order_var);

do_replace=@replacement; %#ok<NASGU>

express=['\<',parser.cell2matize(obj.endogenous.name),'\>'];

for iprob=1:size(inner_probs,1)
    
    prob=regexprep(inner_probs{iprob,2},express,'${do_replace($1)}');
    
    [~,markovChainLocs]=decompose_name(inner_probs{iprob,1});
    
    disp(prob)
    
    inner_probs{iprob,2}=struct('func',str2func(['@(x)',prob]),...
        'chain_loc',markovChainLocs);
    
end

    function out=replacement(vname)
        
        ploc=int2str(find(strcmp(vname,ordered_names)));
        
        out=['x(',ploc,',:)'];
        
    end

    function [currState,markovChainLocs]=decompose_name(tp_name)
        
        underscore=find(tp_name=='_');
        
        currState=str2double(tp_name(underscore(2)+1:underscore(3)-1));
        
        chain_name=tp_name(1:underscore(1)-1);
        
        chain_loc= strcmp(chain_name,obj.markov_chains.chain_names);
        
        markovChainLocs = regimes(:,chain_loc) == currState;
        
    end

end