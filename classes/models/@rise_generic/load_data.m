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


function [obj,issue,retcode]=load_data_intern(obj,varargin)
issue='';
retcode=0;
obj=set(obj,varargin{:});
% check wether the data are provided
obj.options.data=ts.collect(obj.options.data);
data_provided=obj.options.data.NumberOfVariables>0;

if data_provided
    % the length/number of pages of the dataset depends on the horizon of the
    % shocks
    pages=[]; % get all the pages
    [verdier,obj.options.estim_start_date,obj.options.estim_end_date]=...
        utils.time_series.data_request(obj.options.data,obj.observables.name,...
        obj.options.estim_start_date,obj.options.estim_end_date,pages);
    
    % load exogenous and endogenous
    %------------------------------
    is_endogenous=obj.observables.is_endogenous;
    obj.data.y=verdier(is_endogenous,:,:);
    obj.data.x=verdier(~is_endogenous,:,:);
    if obj.options.data_demean
        obj.data.y=bsxfun(@minus,obj.data.y,utils.stat.nanmean(obj.data.y(:,:,1),2));
    end
    obj.data.nobs=size(verdier,2);
    obj.data.npages=size(verdier,3);
    obj.data.start=1;
    obj.data.finish=obj.data.nobs;
    obj.data.varobs_id=real(obj.observables.state_id(is_endogenous));
    
    % conditional forecasting
    %-------------------------
    % Is it restrictive to impose that the variables be observed?
    estim_cond_vars=obj.options.estim_cond_vars;
    if ~isempty(estim_cond_vars)
        if ischar(estim_cond_vars)
            estim_cond_vars=cellstr(estim_cond_vars);
        end
        pos=locate_variables(estim_cond_vars,obj.observables.name,true);
        if any(isnan(pos))
            disp(estim_cond_vars(isnan(pos)))
            error('the variables above are not declared as observables')
        end
        endo_pos=is_endogenous(pos);
        % separate exogenous from exogenous
        %-----------------------------------
        obj.data.restr_y_id=obj.observables.state_id(pos(endo_pos));
        obj.data.restr_x_id=obj.observables.state_id(pos(~endo_pos));
    end
    % add further description fields
    %-------------------------------
    obj.data=data_description(obj.data);
    obj.data_are_loaded=true;
    
else
    retcode=500;
    disp([mfilename,'(GENTLE WARNING):: no actual or simulated data provided for filtering/estimation'])
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


