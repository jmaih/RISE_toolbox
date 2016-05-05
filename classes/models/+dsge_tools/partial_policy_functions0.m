function polfun=partial_policy_functions(m,histdb,state_ranges)

[m,retcode]=solve(m);

if retcode
    
    error(['model could not be solved: ',decipher(retcode)])
    
end

state_list=fieldnames(state_ranges);

polfun=struct();

for ivar=1:numel(state_list)
    
    name=state_list{ivar};
    
    irange=state_ranges.(name);
    
    npoints=numel(irange);
    
    db=histdb;
    
    endhist=db.(name).start;
    
    for ipoint=1:npoints
        
        db.(name)(endhist)=irange(ipoint);
        
        sdev=simulate(m,'simul_historical_data',db,...
            'simul_periods',1,...
            'simul_shock_uncertainty',false,...
            'simul_to_time_series',false);
        
        if ipoint==1
            
            tank=sdev{1}(2*ones(1,npoints),:);
            
            endo_names=sdev{2}{2,2};
            
            nv=numel(endo_names);
            
        end
        
        tank(ipoint,:)=sdev{1}(2,:);
        
        fprintf(1,'state (%s):: point %0.0f of %0.0f\n',name,ipoint,npoints);
    end
    
    tank=mat2cell(tank.',ones(nv,1),npoints);
    
    polfun.(name)=cell2struct(tank,endo_names,1);
    
    polfun.(name).x_axis=irange;
    
end