function m=parameterize(m,fcn)

if isempty(m)
    
    m=cell(0,4);
    
    return
    
end

if nargin<2,fcn=[]; end

if isempty(fcn),fcn=@randn; end

if ischar(fcn),fcn=str2func(fcn);end

m.parameter_values=set_parameters();

m=set_auxiliary_switching_parameters(m);


    function p=set_parameters()
        
        p=m.parameter_values;
        
        plist=get(m,'par_list');
        
        p=fcn(size(p));
        
        chains=m.markov_chains.chain_names-'const';
        
        if isempty(chains)
            
            return
            
        end
        
        regs=m.markov_chains.regimes(:,2:end);
        
        for ichain=1:numel(chains)
            
            cname=chains{ichain};
            
            loc=strcmp(regs(1,:),cname);
            
            states=unique(cell2mat(regs(2:end,loc)));
            
            nstates=numel(states);
            
            for s0=1:nstates
                
                p01=rand(1,nstates);
                
                p01=p01/sum(p01);
                
                for s1=1:nstates
                    
                    pname=sprintf('%s_tp_%0.0f_%0.0f',cname,s0,s1);
                    
                    ploc=find(strcmp(plist,pname));
                    
                    if isempty(ploc),continue,end
                    
                    p(ploc,:)=p01(s1);
                    
                end
                
            end
            
        end
        
    end

end