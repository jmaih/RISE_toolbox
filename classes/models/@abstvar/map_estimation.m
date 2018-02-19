function [estimList,start,elb,eub,theMap,transProbs]=...
    map_estimation(markov_chains,regimes,p)

    function f=map_one_parameter(ip,loc,states)
        
        ns=numel(states);
        
        f=[ip*ones(ns,1),...
            loc*ones(ns,1),states(:)];
        
    end

nparams=numel(p.pnames);

subregs=cell2mat(regimes(2:end,2:end));

nregimes=size(subregs,1);

N=10000;

estimList=cell(N,1); start=nan(N,1); elb=nan(N,1); eub=nan(N,1);

theMap=nan(N,3);

iter=0;

iter_map=0;

ip=0;

for ic=1:numel(markov_chains)
    
    nstates=markov_chains(ic).number_of_states;
    
    ch_name=markov_chains(ic).name;
    
    param_list=markov_chains(ic).param_list;
    
    pos=locate_variables(param_list,p.pnames);
    
    if nstates==1
        
        do_constant_pocket(param_list,pos)
        
    else
        
        do_one_switching_pocket(ch_name,param_list,pos)
        
    end
    
end

estimList=estimList(1:iter);

start=start(1:iter);

elb=elb(1:iter);

eub=eub(1:iter);

% estimParams,rowsInParams,colsInParams
theMap=abstvar.encode_map(theMap(1:iter_map,:),[nparams,nregimes]);

transProbs=do_parameter_transform();

    function do_one_switching_pocket(ch_name,param_list,pos)
        
        V=cell(1,nstates);
        
        tmp_regimes=V;
        
        for irow=1:numel(param_list)
            
            param=param_list{irow};
            
            pos_irow=pos(irow);
            
            for istate=1:nstates
                
                if irow==1
                    
                    tmp_regimes{istate}=find(subregs(:,ic)==istate);
                    
                end
                
                V{istate}=sprintf('%s_%s_%0.0f',param,ch_name,istate);
                
                ip=ip+1; % estimated parameter index
                
                theseRegimes=tmp_regimes{istate};
                
                nsubregimes=numel(theseRegimes);
                
                theMap(iter_map+(1:nsubregimes),:)=map_one_parameter(ip,pos_irow,...
                    theseRegimes); % ip, pos(ip),1:nregimes
                
                iter_map=iter_map+nsubregimes;
                
            end
            
            rr=iter+(1:nstates);
            
            estimList(rr)=V(:);
            
            elb(rr)=p.lb(pos_irow)*ones(nstates,1);
            
            eub(rr)=p.ub(pos_irow)*ones(nstates,1);
            
            start(rr)=p.start(pos_irow)*ones(nstates,1);
            
            iter=rr(end);
            
        end
        
    end

    function do_constant_pocket(param_list,pos)
        % constant
        
        np1=numel(param_list);
        
        prange=iter+(1:np1);
        
        estimList(prange)=param_list(:);
        
        elb(prange)=p.lb(pos);
        
        eub(prange)=p.ub(pos);
        
        start(prange)=p.start(pos);
        
        for ipos=1:np1
            
            ip=ip+1; % estimated parameter index
            
            theMap(iter_map+(1:nregimes),:)=map_one_parameter(ip,pos(ipos),...
                (1:nregimes)); % ip, pos(ip),1:nregimes
            
            iter_map=iter_map+nregimes;
            
        end
        
        iter=prange(end);
        
    end

    function transProbs=do_parameter_transform()
        
        % transform each dirichlet
        transProbs=cell(1,1000);
        
        iter_prob=0;
        
        for ichain=2:numel(markov_chains)
            
            chain=markov_chains(ichain);
            
            cname=chain.name;
            
            n=chain.number_of_states;
            
            if n<=2
                % no need to transform single probabilities
                continue
                
            end
            
            for ii=1:n
                
                stud=sprintf('%s_tp_%0.0f_\\d+',cname,ii);
                
                collection=regexp(estimList,stud,'match');
                
                collection=[collection{:}];
                
                iter_prob=iter_prob+1;
                
                transProbs{iter_prob}=locate_variables(collection,estimList);
                
                % what to do about the covariances? They need to be
                % positive definite and inverse-wisharing them here is
                % not much of a help since some of their elements may
                % switch.
                
            end
            
        end
        
        transProbs(iter_prob+1:end)=[];
        
    end

end

