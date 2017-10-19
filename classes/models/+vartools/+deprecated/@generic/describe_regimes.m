function the_regimes=describe_regimes(markov_chain_info)

the_regimes=markov_chain_info.regimes;

nregs=markov_chain_info.regimes_number;

nchains=markov_chain_info.chains_number;

chain_names=markov_chain_info.chain_names;

tmp=cell(1,nregs);

for ireg=1:nregs
    
    for ichain=1:nchains
        
        new_item=[chain_names{ichain},' = ',sprintf('%0.0f',the_regimes{ireg+1,ichain+1})];
        
        if ichain==1
            
            tmp{ireg}=new_item;
            
        else
            
            tmp{ireg}=[tmp{ireg},' & ',new_item];
            
        end
        
    end
    
end

the_regimes=tmp;

end
