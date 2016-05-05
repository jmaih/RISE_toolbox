function polfun=partial_policy_functions(m,histdb,state_ranges,istate)

if nargin<4
    
    istate=[];
    
    if nargin<3
        
        state_ranges=[];
        
        if nargin <2
            
            histdb=[];
            
        end
        
    end
    
end

[m,retcode]=solve(m);

if retcode
    
    error(['model could not be solved: ',decipher(retcode)])
    
end

endo_list=get(m,'endo_list');

nv=numel(endo_list);

[y0_iov,T,sstate,new_order,iov,xloc,log_vars,sig,shocks]=initialize();

state_list=fieldnames(state_ranges);

polfun=struct();

for ivar=1:numel(state_list)
    
    name=state_list{ivar};
    
    name_pos=strcmp(name,endo_list);
    
    y0i=y0_iov;
    
    irange=state_ranges.(name);
    
    npoints=numel(irange);
    
    for ipoint=1:npoints
        
        y0i(name_pos)=irange(ipoint);
        
        y1=simulate_one_step();
        
        if ipoint==1
            
            tank=y1(:,ones(1,npoints));
            
        end
        
        tank(:,ipoint)=y1;
        
        fprintf(1,'state (%s):: point %0.0f of %0.0f\n',name,ipoint,npoints);
    end
    
    tank=mat2cell(tank,ones(nv,1),npoints);
    
    polfun.(name)=cell2struct(tank,endo_list,1);
    
    polfun.(name).x_axis=irange;
end

    function [yiov,T,sstate,new_order,iov,xloc,log_vars,sig,shocks]=initialize()
        
        [T,~,sstate,new_order,xloc,log_vars]=load_solution(m,'ov',false);
        
        sig=1;
        
        iov(new_order)=1:numel(new_order);
        
        if isempty(istate)
            
            istate=1;
            
        end
        
        if isempty(histdb)
            
            yiov=sstate{istate}(iov);
            
        else
            
            yiov=nan(size(sstate{1}));
            
            for ilist=1:numel(endo_list)
                
                yiov(ilist)=double(histdb.(endo_list{ilist}));
                
            end
            
        end
        
        nshocks=sum(m.exogenous.number);
        
        horizon=1+max(m.exogenous.shock_horizon);
        
        shocks=zeros(nshocks,horizon);
        
        if isempty(state_ranges)
            
            state_list=initialize_state_ranges(m);
            
        end
        
    end

    function y1=simulate_one_step()
        
        y00=y0i(new_order);
        
        y00(log_vars)=log(y00(log_vars));
        
        y0=struct('y',y00);
        
        y1=utils.forecast.one_step_engine(T(:,istate),y0,sstate{:,istate},xloc,sig,shocks);
        
        y1.y(log_vars)=exp(y1.y(log_vars));
        
        y1=y1.y(iov);
        
    end

end

function state_list=initialize_state_ranges(m)
% verify the list of state variables and get sstate
state_list=get(m,'endo_list(state)');

% get steady state
ss=get(m,'sstate');

% define range for state variables
delta=0.3;

npoints=300;

state_range=struct();

for iv=1:numel(state_list)
    
    name=state_list{iv};
    
    ssi=ss.(name);
    
    state_range.(name)=linspace(ssi*delta,ssi*1/delta,npoints);
end

end