function [obj,issue,retcode]=load_data(obj,varargin)%,estimation_flag
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


if isempty(obj)
    obj=struct('data',ts,...
        'data_demean',false);
    % data for conditional forecast and real-time estimation
    return
end
nobj=numel(obj);
issue='';
retcode=0;
for ii=1:nobj
    if isa(obj(ii),'dsge') && obj(ii).is_optimal_simple_rule_model
        if obj(ii).observables.number(1)>0
            warning([mfilename,':: data not needed OPTIMAL SIMPLE RULES problems'])
        end
        obj(ii).data_are_loaded=true;
        continue
    else
        if obj(ii).observables.number(1)==0
            error([mfilename,':: List of observables not provided for model ''',obj(ii).filename,''])
        end
        [obj(ii),issue,rcode]=load_data_intern(obj(ii),varargin{:});
        if nobj>1 && ~isempty(issue)
            warning([mfilename,':: problem for model ''',obj(ii).filename,''': ',issue]) %#ok<WNTAG>
        end
        retcode=max(retcode,rcode);
    end
end

end


function [obj,issue,retcode]=load_data_intern(obj,varargin)
issue='';
retcode=0;
obj=set(obj,varargin{:});
% check wether the data are provided
obj.options.data=ts.collect(obj.options.data);
data_provided=obj.options.data.NumberOfVariables>0;

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
    
    % the length/number of pages of the dataset depends on the horizon of the
    % shocks but for the moment, we consider all pages
    pages=[];
    [verdier,obj.options.estim_start_date,obj.options.estim_end_date]=...
        utils.time_series.data_request(obj.options.data,allvars,...
        obj.options.estim_start_date,obj.options.estim_end_date,pages);
    
    obs_id=locate_variables(obj.observables.name,allvars);
    
    % load endogenous and deterministic exogenous
    %----------------------------------------------
    obj.data.y=verdier(obs_id(is_endogenous_obs),:,:);
    obj.data.x=verdier(obs_id(~is_endogenous_obs),:,:);
    if obj.options.data_demean
        obj.data.y=bsxfun(@minus,obj.data.y,utils.stat.nanmean(obj.data.y(:,:,1),2));
    end
    obj.data.nobs=size(verdier,2);
    obj.data.npages=size(verdier,3);
    obj.data.start=1;
    obj.data.finish=obj.data.nobs;
    obj.data.varobs_id=real(obj.observables.state_id(is_endogenous_obs));
    
    % conditional forecasting data
    %------------------------------
    Locb = locate_variables(forecast_cond_exo_vars,allvars,true);
    obj.data.z=verdier(Locb,:,:);
    
    if isempty(pos_endo_fkst_in_obs)
        obj.data.restr_y_id=[];
    else
        % id in state vector + (id in the observables)*1i 
        %--------------------------------------------------
        obj.data.restr_y_id=obj.observables.state_id(pos_endo_fkst_in_obs)+...
            pos_endo_fkst_in_obs*1i;
    end
    
    if isempty(pos_exo_fkst_in_exo)
        obj.data.restr_z_id=[];
    else
        obj.data.restr_z_id=pos_exo_fkst_in_exo+pos_exo_fkst_in_exo*1i;
    end
    
    % add further description fields
    %-------------------------------
    obj.data=data_description(obj.data);
    obj.data_are_loaded=true;
else
    retcode=500;
    disp([mfilename,'(GENTLE WARNING):: no actual or simulated data provided for filtering/estimation'])
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

% points to include in the calculation of the likelihood
%-------------------------------------------------------
d.include_in_likelihood=false(1,smpl);
d.include_in_likelihood(d.start:d.finish)=true;
end


