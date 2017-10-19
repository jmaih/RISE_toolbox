function [dataOut,dataset,estim_start_date,estim_end_date,retcode]=data_prerequest(obj,dataset)

% check wether the data are provided
dataset=ts.collect(dataset);

data_provided=dataset.NumberOfVariables>0;

dataOut=[];

retcode=0;
    
%     obj.data_are_loaded=true;

if data_provided
    
    is_endogenous_obs=obj.observables.is_endogenous;
    
    xvars=obj.observables.name(~is_endogenous_obs);
    
    % list of all the variables: union(observables,cond_endo,cond_exo)
    %------------------------------------------------------------------
    [pos_endo_fkst_in_obs,forecast_cond_endo_vars]=find_positions(obj.options.forecast_cond_endo_vars,'observables');
    % exclude the deterministic exogenous from the list and thereby ensure
    % there is no collision between the deterministic exogenous and the
    % stochastic exogenous 
    [pos_exo_fkst_in_exo,forecast_cond_exo_vars]=find_positions(obj.options.forecast_cond_exo_vars,'exogenous',xvars);

    allvars=union(obj.observables.name,...
        union(forecast_cond_endo_vars,forecast_cond_exo_vars));
    
    % remove empty cells to avoid a warning about not-found variables
    %------------------------------------------------------------------
    allvars=allvars(cellfun(@(x)~isempty(x),allvars,'uniformOutput',true));
    
    % the length/number of pages of the dataset depends on the horizon of the
    % shocks but for the moment, we consider all pages
    pages=[];
    
    [verdier,estim_start_date,estim_end_date]=...
        utils.time_series.data_request(dataset,allvars,...
        obj.options.estim_start_date,obj.options.estim_end_date,pages);
    
    % locate the log vars and log them
    %----------------------------------
    if isa(obj,'dsge')
        
        verdier=do_log(verdier,allvars);
        
    end
    
    
    obs_id=locate_variables(obj.observables.name,allvars);
    
    % load endogenous and deterministic exogenous
    %----------------------------------------------
    dataOut=struct();
    
    dataOut.y=verdier(obs_id(is_endogenous_obs),:,:);
    
    dataOut.x=verdier(obs_id(~is_endogenous_obs),:,:);

	dataOut.ymean=[];
    
    if obj.options.data_demean
        
		dataOut.ymean=utils.stat.nanmean(dataOut.y(:,:,1),2);
    
        dataOut.y=bsxfun(@minus,dataOut.y,dataOut.ymean);
        
    end
    
    dataOut.nobs=size(verdier,2);
    
    dataOut.npages=size(verdier,3);
    
    dataOut.start=1;
    
    dataOut.finish=dataOut.nobs;
    
    dataOut.varobs_id=real(obj.observables.state_id(is_endogenous_obs));
    
    % conditional forecasting data
    %------------------------------
    Locb = locate_variables(forecast_cond_exo_vars,allvars,true);
    
    dataOut.z=verdier(Locb,:,:);
    
    if isempty(pos_endo_fkst_in_obs)
        
        dataOut.restr_y_id=[];
        
    else
        
        % id in state vector + (id in the observables)*1i 
        %--------------------------------------------------
        dataOut.restr_y_id=obj.observables.state_id(pos_endo_fkst_in_obs)+...
            pos_endo_fkst_in_obs*1i;
        
    end
    
    if isempty(pos_exo_fkst_in_exo)
        
        dataOut.restr_z_id=[];
        
    else
        
        dataOut.restr_z_id=pos_exo_fkst_in_exo+pos_exo_fkst_in_exo*1i;
        
    end
    
    % add further description fields
    %-------------------------------
    dataOut=data_description(dataOut);
    
else
    
    retcode=500;
    
    disp([mfilename,'(GENTLE WARNING):: no actual or simulated data provided for filtering/estimation'])
    
end


    function v=do_log(v,vnames)
        
        is_log_var=obj.log_vars;
        
        silent=true;
        
        locs=locate_variables(vnames,obj.endogenous.name,silent);
        
        for ivar=1:numel(vnames)
            
            if isnan(locs(ivar))
                
                continue
                
            end
            
            if is_log_var(locs(ivar))
                
                v(ivar,:,:)=log(v(ivar,:,:));
                
            end
        end
    end


    function [pos,condvars]=find_positions(condvars,type,exclude)
        
        if nargin<3
            
            exclude={};
            
        end
        
        pos=[];
        
        if ~isempty(condvars)
            
            if ischar(condvars)
                
                condvars=cellstr(condvars);
                
            end
            
            main_names=setdiff(obj.(type).name,exclude);
            
            pos=locate_variables(condvars,main_names,true);
            
            if any(isnan(pos))
                
                disp(condvars(isnan(pos)))
                
                if strcmp(type,'exogenous')
                    
                    disp(['Note that deterministic exogenous variables',...
                        ' cannot be declared in the forecast_cond_exo_vars group'])
                end
                
                error(['the conditional variables above are not declared as ',type])
                
            end
            
        end
        
    end

end

function [d]=data_description(d)

smpl=size(d.y,2);

% expand the pages into one page
d.data_structure=permute(d.y,[2,1,3]);

d.data_structure=permute(d.data_structure(:,:),[2,1]);

d.data_structure=~isnan(d.data_structure);

tmp=find(any(d.data_structure(:,1:d.finish,1)==false,1),1,'last');

if isempty(tmp)
    
    tmp=0;
    
end

d.no_more_missing=min(d.finish,tmp+1);

% re-fold
d.data_structure=~isnan(d.y);

% conditional data
%------------------
dz=permute(d.z,[1,3,2]);

last_non_nan=zeros(1,smpl);

for t=1:smpl
    
    good=find(any(~isnan(dz(:,:,t)),1),1,'last');
    
    if isempty(good)
        
        good=0;
        
    end
    
    last_non_nan(t)=good;
    
end

d.last_good_conditional_observation=last_non_nan;

end