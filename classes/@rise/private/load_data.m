function [obj,issue,retcode]=load_data(obj,varargin)%,estimation_flag

% if nargin<2
%     estimation_flag=true;
% end
nobj=numel(obj);
issue='';
retcode=0;
for ii=1:nobj
    if obj(ii).is_optimal_simple_rule_model 
        if obj(ii).NumberOfObservables(1)>0
            error([mfilename,':: Observables are not allowed with OPTIMAL SIMPLE RULES problems'])
        end
    else
        if obj(ii).NumberOfObservables(1)==0
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
obj.options=mysetfield(obj.options,varargin{:});
% check wether the data are provided
data_provided=obj.options.data.NumberOfVariables>0;
% check whether there is a simulation
simvals=vertcat(obj.varendo.value);
simulation_available=~isempty(simvals);
% if ~data_provided && ~simulation_available
%     if estimation_flag
%         error([mfilename,':: no data available and no simulation provided as replacement'])
%     else
%         issue='data not provided';
%         return
%     end
% end

% now we can load the data
% check whether there are conditional data
% collect all the conditioning variables
cond_names0=cell(0);
if obj.options.cond_data_ct.NumberOfObservations
    cond_names0=union(cond_names0,cellstr(obj.options.cond_data_ct.varnames));
end
if obj.options.cond_data_lb.NumberOfObservations
    cond_names0=union(cond_names0,cellstr(obj.options.cond_data_lb.varnames));
end
if obj.options.cond_data_ub.NumberOfObservations
    cond_names0=union(cond_names0,cellstr(obj.options.cond_data_ub.varnames));
end
ncv0=numel(cond_names0);
%     options.cond_data_ct=rise_time_series;
%     options.cond_data_lb=rise_time_series;
%     options.cond_data_ub=rise_time_series;

if data_provided
    The_data=rise_time_series.collect(obj.options.data);
    if isempty(obj.options.estim_start_date)
        obj.options.estim_start_date=The_data.start;
    end
    if isempty(obj.options.estim_end_date)
        obj.options.estim_end_date=The_data.finish;
    end
    start=The_data.TimeInfo(1).date_2_observation(obj.options.estim_start_date);
    finish=The_data.TimeInfo(1).date_2_observation(obj.options.estim_end_date);
    if start<=0
        error('estimation start date inferior to data start')
    end
    if finish>The_data.NumberOfObservations
        error('estimation end date greater than to data end date')
    end
    smpl=finish-start+1;
    ids=locate_variables({obj.varobs.name},The_data.varnames);
    DataValues=double(The_data);
    verdier=transpose(DataValues(start:finish,ids));
    if obj.options.data_demean
        verdier=bsxfun(@minus,verdier,mean(verdier,2));
    end
    verdier=mat2cell(verdier,ones(obj.NumberOfObservables(1),1),smpl);

    % load all the data if possible, in case we do not want to
    % update the state for each parameter when doing forecasting
    ids_all=locate_variables({obj.varendo.name},The_data.varnames,true);
    if ~any(isnan(ids_all))
        tmp=transpose(DataValues(start:finish,ids_all));
        for ii=1:obj.NumberOfEndogenous(2)
            obj.varendo(ii)=obj.varendo(ii).set_properties('value',tmp(ii,:));
        end
        clear tmp
    end

    % load exogenous observables
    if obj.NumberOfObservables(2)
        ids_all=locate_variables({obj.varobs_exo.name},The_data.varnames,false);
        % scream if you don't find a variable
        tmp=transpose(DataValues(start:finish,ids_all));
        % scream even louder if there is a nan in the exogenous coz then
        % estimation cannot take place
% % %         if any(any(isnan(tmp)))
% % %             error([mfilename,':: missing values in the exogenous, estimation cannot take place'])
% % %         end
        for ii=1:obj.NumberOfObservables(2)
            obj.varobs_exo(ii)=obj.varobs_exo(ii).set_properties('value',tmp(ii,:));
        end
        clear tmp
    end
    
    if ncv0
        [cond_data_ct,horiz_ct,date_start_ct]=dispatch_conditional_data(obj.options.cond_data_ct);
        [cond_data_lb,horiz_lb,date_start_lb]=dispatch_conditional_data(obj.options.cond_data_lb);
        [cond_data_ub,horiz_ub,date_start_ub]=dispatch_conditional_data(obj.options.cond_data_ub);
        
        horiz0=max([horiz_ct,horiz_lb,horiz_ub]);
        % the expansion order should be less than or equal to the horizon of the conditional
        % data
        if obj.options.order>horiz0
            error([mfilename,':: expansion order cannot be greater than the horizon of the conditioning data'])
        end
        
        % now get the conditional data for the the central tendency
        cond_data_ct=format_conditional_data(cond_data_ct,...
            date_start_ct,cond_names0,start,horiz0);
        
        % do the same with the lower bound and the upper bound
        cond_data_lb=format_conditional_data(cond_data_lb,...
            date_start_lb,cond_names0,start,horiz0);
        
        cond_data_ub=format_conditional_data(cond_data_ub,...
            date_start_ub,cond_names0,start,horiz0);
        
        % Now that all these databases start at the same date, we can
        % harmonize their span
        span_ct=size(cond_data_ct,3);
        span_lb=size(cond_data_lb,3);
        span_ub=size(cond_data_ub,3);
        span0=max(span_ct,max(span_lb,span_ub));
        cond_data_ct=cat(3,cond_data_ct,nan(ncv0,horiz0,span0-span_ct));
        cond_data_lb=cat(3,cond_data_lb,nan(ncv0,horiz0,span0-span_lb));
        cond_data_ub=cat(3,cond_data_ub,nan(ncv0,horiz0,span0-span_ub));
        % now check that the bounds are consistent
        tmp=cond_data_ct-cond_data_lb;
        good=~isnan(tmp);
        if any(tmp(good)<0)
            error([mfilename,':: some central-tendency elements lower than their lower-bound counterpart'])
        end
        tmp=cond_data_ub-cond_data_ct;
        good=~isnan(tmp);
        if any(tmp(good)<0)
            error([mfilename,':: some central-tendency elements greater than their upper-bound counterpart'])
        end
        
        % now construct the cond_varobs object
       obj.NumberOfConditionalObservables=ncv0;
       locs=locate_variables(cond_names0,{obj.varendo.name});
       for ii=1:ncv0
           tex_name=obj.varendo(locs(ii)).tex_name;
           obj.cond_varobs(ii,1)=rise_variable(cond_names0{ii},...
               'tex_name',tex_name,'id',locs(ii),...
               'value',cond_data_ct(ii,:,:),...
               'lb',cond_data_lb(ii,:,:),'ub',cond_data_ub(ii,:,:));
       end
        % REMAINS: CONDITIONING ON SHOCKS
    end
elseif simulation_available
    ids=locate_variables({obj.varobs.name},{obj.varendo.name});
    verdier=mat2cell(simvals(ids,:),...
        ones(1,obj.NumberOfObservables(1)),size(simvals,2));
    if ncv0 % && estimation_flag
        error([mfilename,':: please provide some actual data to estimate with real-time information'])
    end
else
    retcode=500;
    disp([mfilename,'(GENTLE WARNING):: no actual or simulated data provided for filtering/estimation'])
end

if simulation_available || data_provided
%     [obj.varobs(:).value]=(verdier{:}); % this vectorized form works but
%     does not update the other fields in the objects
    for ii=1:obj.NumberOfObservables
        obj.varobs(ii)=obj.varobs(ii).set_properties('value',verdier{ii});
    end
    if simulation_available && ~data_provided %&& estimation_flag 
        issue='simulated data used for estimation in lieu of actual data. Starting date will be 2000Q1';
        warning([mfilename,':: ',issue]) %#ok<WNTAG>
		estim_start_date='2000Q1';
	else
		estim_start_date=obj.options.estim_start_date;
    end
	nobs=size(obj.varobs(1).value,2);
    obj.dates_filtering=rise_date(estim_start_date)+(0:nobs);
    obj.dates_smoothing=obj.dates_filtering(1:end-1);

    if obj.is_dsge_var_model
        data=vertcat(obj.varobs.value)';
        const=obj.options.dsge_var_constant;
        n=obj.NumberOfObservables;
        p=obj.options.dsge_varlag;
        Y=data(p+1:end,:);
        smpl=size(Y,1);
        
        X=nan(smpl,const+n*p);
        if const
            X(:,1)=1;
        end
        for ii=1:p
            X(:,const+((ii-1)*n+1:ii*n))=data(p+1-ii:end-ii,:);
        end
        dsge_var=struct('YY',Y'*Y,'YX',Y'*X,...
            'XX',X'*X,'XY',X'*Y,'T',smpl,...
            'n',n,'p',p,'constant',const);
        obj=obj.set_properties('dsge_var',dsge_var);
    end
    obj.data_are_loaded=true;
 
    % Provision for Endogenous Priors
    if ~isempty(obj.options.prior_endogenous)
        prior_endogenous=obj.options.prior_endogenous;
        targets=[];
        prior_sample=[];
        do_it=false;
        if islogical(prior_endogenous) 
            if prior_endogenous==true
                do_it=true;
            else
                obj.endogenous_priors=[];
            end
        elseif isstruct(prior_endogenous)
            if isfield(prior_endogenous,'targets')
                targets=prior_endogenous.targets;
                do_it=true;
            end
            if isfield(prior_endogenous,'prior_sample')
                prior_sample=prior_endogenous.prior_sample;
                do_it=true;
            end
        end
        if do_it
            obj.endogenous_priors=rise_endo_priors(obj,targets,prior_sample);
        end
    end
elseif obj.is_optimal_simple_rule_model 
    obj.data_are_loaded=true;
end

function [cond_data,horiz,date_start]=dispatch_conditional_data(cond_data)
horiz=0;
date_start=[];
if ~isempty(cond_data)
    cond_data=rise_time_series.collect(cond_data);
    horiz=cond_data.NumberOfPages;
    date_start=cond_data.TimeInfo(1);
else
    cond_data=[];
end



function cond_data=format_conditional_data(cond_data_ts,The_data_TimeInfo,cond_names0,start,horiz0)
ncv0=numel(cond_names0);
if isempty(cond_data_ts)
    cond_data=nan(ncv0,horiz0,1);
else
    cond_names=cellstr(cond_data_ts.varnames);
    ncv=cond_data_ts.NumberOfVariables;
    locs=locate_variables(cond_names,cond_names0);
    horiz=cond_data_ts.NumberOfPages;
    % if ~all(ismember(cond_names,cond_names0))
    %     error([mfilename,':: all conditioning databases should have the same variables'])
    % end
    if cond_data_ts.TimeInfo(end)<The_data_TimeInfo(1)
        error([mfilename,':: conditional information cannot end before the start of the actual data'])
    end
    add_on_beg=0;add_on_end=0;discard_beg=0;
    if cond_data_ts.TimeInfo(1)<The_data_TimeInfo(1)
        % cut out the difference at the beginning
        discard_beg=numel(cond_data_ts.TimeInfo(1):The_data_TimeInfo(1)-1);
    elseif cond_data_ts.TimeInfo(1)>The_data_TimeInfo(1)
        % fill in the missing at the beginning
        add_on_beg=numel(The_data_TimeInfo(1):cond_data_ts.TimeInfo(1)-1);
    end
    cond_data=cell2mat(cond_data_ts.data(2:end,2:end,:));
    % permute to names x horizon x smpl
    cond_data=permute(cond_data,[2,3,1]);
    cond_data=cond_data(:,:,discard_beg+1:end);
    
    cond_data=cat(3,nan(ncv,horiz,add_on_beg),cond_data,nan(ncv,horiz,add_on_end));
    % now we trim the beginning to match the estimation data
    if start<size(cond_data,3)
        cond_data=cond_data(:,:,start:end);
    else
        error([mfilename,':: conditional information cannot end before the start of the estimation sample'])
    end
    % We do not cut the end since the information is potentially useful
    % for forecasting
    % now put everything into the biggest format
    tmp=cond_data;
    cond_data=nan(ncv0,horiz0,size(tmp,3));
    cond_data(locs,1:horiz,:)=tmp;
end
