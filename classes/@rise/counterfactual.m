function obj=counterfactual(obj,param_type,shocks_db,param_struct)
% shocks_db is a rise_time_series object with alternative history for the shocks
if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    obj=struct();
    return
end
if nargin < 4
    param_struct=[];
    if nargin<3
        shocks_db=[];
        if nargin<2
            param_type=[];
        end
    end
end

mode_=vertcat(obj.estimated_parameters.mode);
mean_=vertcat(obj.estimated_parameters.mean);
startval=vertcat(obj.estimated_parameters.startval);
params=[];
if isempty(param_type)
    if ~isempty(mean_)
        params=mean_;
    elseif ~isempty(mode_)
        params=mode_;
    end
else
    if strcmp(param_type,'mode')
        if isempty(mode_)
            warning([mfilename,':: mode has not been estimated']) %#ok<WNTAG>
        end
        params=mode_;
    elseif strcmp(param_type,'mean')
        if isempty(mean_)
            if isempty(mode_)
                warning([mfilename,':: model has not been estimated, using calibration']) %#ok<WNTAG>
            else
                warning([mfilename,':: mean has not been estimated, using the mode']) %#ok<WNTAG>
            end
            params=mode_;
        else
            params=mean_;
        end
    elseif strcmp(param_type,'calibration')
        params=startval;
    else
        error([mfilename,':: unrecognized parameter type ',param_type])
    end
end

% compute the history
[obj,~,~,retcode]=filter(obj,'evaluate_params',params,'kf_filtering_level',3);
if retcode
    error([mfilename,':: model could not be solved for this parameterization'])
end

if ~isempty(shocks_db)
    T=obj.T;
    R=obj.R; %
    [endo_nbr,~,order,regime_nbr]=size(R);
    if regime_nbr>1
        error([mfilename,':: counterfactual analysis for multiple-regime models not yet implemented'])
    end
    SS=obj.steady_state_and_balanced_growth_path(1:endo_nbr,:);
    % use the new shocks and the smoothed information to reconstruct history
    
    % dates for beginning and end of history
    start_date=obj.Filters.smoothed_shocks.regime_1.TimeInfo(1);
    finish_date=obj.Filters.smoothed_shocks.regime_1.TimeInfo(end);
    % check dates consistency
    if start_date<shocks_db.TimeInfo(1)
        error([mfilename,':: alternative shock history cannot start before ',start_date.date])
    end
    if finish_date<shocks_db.TimeInfo(end)
        error([mfilename,':: alternative shock history cannot end after ',finish_date.date])
    end
    if obj.Filters.smoothed_shocks.regime_1.NumberOfPages<shocks_db.NumberOfPages
        error([mfilename,':: alternative shock history cannot have more pages than the baseline'])
    end
    % locations of the alternative history
    alt_start=start_date.date_2_observation(shocks_db.TimeInfo(1));
    alt_finish=start_date.date_2_observation(shocks_db.TimeInfo(end));
    
    % locate the new shocks in the baseline database
    name_id=locate_variables(shocks_db.varnames,{obj.varexo.name});

    % get the old shocks (baseline case)
    smoothed_shocks=cell2mat(obj.Filters.smoothed_shocks.regime_1.data(2:end,2:end,:));
    
    % replace with the alternative shocks
    smoothed_shocks(alt_start:alt_finish,name_id,1:shocks_db.NumberOfPages)=...
        cell2mat(shocks_db.data(2:end,2:end,:));
    % now permute the order to names,time,expectations
    smoothed_shocks=permute(smoothed_shocks,[2,1,3]);
    
    % get the number of observations
    NumberOfObservations=obj.Filters.smoothed_shocks.regime_1.NumberOfObservations;
    
    % initialize output
    COUNTER=zeros(endo_nbr,NumberOfObservations);
    % recreate history: things like this could be done in a simulate function
    % but I already have one and this is not the time to mess with it.
    % there is no way to know what y0 is and setting it to be the steady
    % state is only good for filtering, not for smoothing. Hence we impose
    % that the first observation (y1) is given and we compute the rest
    COUNTER(:,1)=cell2mat(obj.Filters.smoothed_variables.regime_1.data(2,2:end));
    % we then start the iterations at 2
    for t=2:NumberOfObservations
        COUNTER(:,t)=SS+T*(COUNTER(:,t-1)-SS);
        for r=1:order
            COUNTER(:,t)=COUNTER(:,t)+R(:,:,r)*smoothed_shocks(:,t,r);
        end
    end
    obj=obj.set_properties('counterfactual_shocks',rise_time_series(start_date.date,transpose(COUNTER),char({obj.varendo.name})));
end

if ~isempty(param_struct)
    % change the parameters and redo the filtering
    % 1- set the parameters
    newobj=obj.set_parameters(param_struct);
    % 2- filter using those parameters but save to a different object
    [newobj,~,~,retcode]=newobj.filter('kf_filtering_level',3);
    if retcode
        error([mfilename,':: model could not be solved for this parameterization'])
    end
    % turn the smoothed to a database
    obj=obj.set_properties('counterfactual_parameters',...
         newobj.Filters.smoothed_variables.regime_1);
end

