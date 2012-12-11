function [counterf,actual]=counterfactual(obj,varargin)
% shocks_db is a rise_time_series object with alternative history for the shocks

% should
Defaults=struct('counterfact_shocks_db','',...
    'counterfact_nullify_others',true);

if isempty(obj)
    counterf=Defaults;
    return
end
nobj=numel(obj);
Fields=fieldnames(Defaults);

counterf=cell(1,nobj);
actual=cell(1,nobj);
for ii=1:nobj
    [counterf{ii},actual{ii}]=counterfactual_intern(obj(ii));
end

if nobj==1
	counterf=counterf{1};
	actual=actual{1};
else
	% pad them as we did with irfs
	
end

    function [COUNTER,ACTUAL]=counterfactual_intern(obj)

        obj=set_options(obj,varargin{:});
        oldoptions=obj.options;
        thisDefault=Defaults;
        for ifield=1:numel(Fields)
            vi=Fields{ifield};
            if isfield(oldoptions,vi)
                thisDefault.(vi)=oldoptions.(vi);
            end
        end
		shocks_db=thisDefault.counterfact_shocks_db;
        counterfact_nullify_others=thisDefault.counterfact_nullify_others;
        
        % get the historical decomposition to partial out the elements that
        % are changed by the counterfactual shocks
		[Histdec,obj]=historical_decomposition(obj);
        
        shockList={obj.varexo.name};
        endoList={obj.varendo.name};
        smoothed_shocks=obj.Filters.Expected_smoothed_shocks;
        smoothed_variables=obj.Filters.Expected_smoothed_variables;
        TimeInfo=smoothed_variables.(endoList{1}).TimeInfo;
        if ischar(shocks_db)
            shocks_db=cellstr(shocks_db);
        end
        if isstruct(shocks_db) || isa(shocks_db,'rise_time_series')
            if isa(shocks_db,'rise_time_series')
                shocks_db=pages2struct(shocks_db);
            end
            if ~isequal(shocks_db.TimeInfo,smoothed_shocks.(shockList{1}).TimeInfo)
                error('shocks database should match history used for filtering')
            end
            % a structure time series
            db_shock_names=fieldnames(shocks_db);
            for ishock=1:numel(shockList)
                shock_name=shockList{ishock};
                if ismember(shock_name,db_shock_names)
                    smoothed_shocks.(shock_name)=shocks_db.(shock_name);
                elseif ismember(shock_name,shockList)
                    smoothed_shocks.(shock_name)=nullify(smoothed_shocks.(shock_name));
                end
            end
        elseif iscellstr(shocks_db)
            % list of shocks to include in the simulation. Get them from
            % the smoother and zero-out all other shocks. but first check
            % that they indeed are present in the list of shocks.
%             locs=locate_variables(shocks_db,shockList);
            other_shocks=setdiff(shockList,shocks_db);
            for ishock=1:numel(other_shocks)
                shock_name=other_shocks{ishock};
                smoothed_shocks.(shock_name)=nullify(smoothed_shocks.(shock_name));
            end
        end
        % get the number of observations
        NumberOfObservations=smoothed_shocks.(shockList{1}).NumberOfObservations;
        
        % now reconstruct the history that would have been
        smoothed_shocks=rise_time_series.collect(smoothed_shocks);
        tmp_names=smoothed_shocks.varnames;
        locs=locate_variables(tmp_names,shockList);
        smoothed_shocks=double(smoothed_shocks);
        smoothed_shocks=permute(smoothed_shocks,[2,1,3]);
        smoothed_shocks=smoothed_shocks(locs,:,:);
        smoothed_variables=rise_time_series.collect(smoothed_variables);
        tmp_names=smoothed_variables.varnames;
        locs=locate_variables(tmp_names,endoList);
        smoothed_variables=double(smoothed_variables);
        smoothed_variables=permute(smoothed_variables,[2,1]);
        smoothed_variables=smoothed_variables(locs,:);
        
        T=obj.T;
        R=obj.R; %
        [endo_nbr,~,order,regime_nbr]=size(R);
        
        SS=obj.steady_state_and_balanced_growth_path(1:endo_nbr,:);
        % use the new shocks and the smoothed information to reconstruct history
        
        % initialize output
        COUNTER=zeros(endo_nbr,NumberOfObservations);
        % recreate history: things like this could be done in a simulate function
        % but I already have one and this is not the time to mess with it.
        % there is no way to know what y0 is and setting it to be the steady
        % state is only good for filtering, not for smoothing. Hence we impose
        % that the first observation (y1) is given and we compute the rest
        COUNTER(:,1)=smoothed_variables(:,1);
        % we then start the iterations at 2
        for t=2:NumberOfObservations
            COUNTER(:,t)=SS+T*(COUNTER(:,t-1)-SS);
            for r=1:order
                COUNTER(:,t)=COUNTER(:,t)+R(:,:,r)*smoothed_shocks(:,t,r);
            end
        end
        
        % put to ...
        COUNTER=double2rise_time_series(COUNTER);
        ACTUAL=double2rise_time_series(smoothed_variables);
        clear smoothed_variables
        
        function dbout=double2rise_time_series(dbin)
            dbout=struct();
            for jj=1:numel(endoList)
                dbout.(endoList{jj})=rise_time_series(TimeInfo,dbin(jj,:)',endoList{jj});
            end
        end
            
        function db=nullify(db)
            if counterfact_nullify_others
                datta=zeros(size(double(db)));
                this_TimeInfo=db.TimeInfo;
                varnames=db.varnames;
                db=rise_time_series(this_TimeInfo,datta,varnames);
            end
        end
    end
end



