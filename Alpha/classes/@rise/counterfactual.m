function [counterf,actual]=counterfactual(obj,varargin)
% shocks_db is a rise_time_series object with alternative history for the shocks

% should
if isempty(obj)
counterf=struct('counterfact_shocks_db','',...
    'counterfact_nullify_others',true,...
    'counterfact_regime',[]); % counterfact_regime={'coef',1','vol',3}
    return
end
nobj=numel(obj);

counterf=cell(1,nobj);
actual=cell(1,nobj);
for ii=1:nobj
    [counterf{ii},actual{ii}]=counterfactual_intern(obj(ii));
end

if nobj==1
	counterf=counterf{1};
	actual=actual{1};
else
	% pad them
	counterf=concatenate_series_from_different_models(counterf);
end

    function [COUNTER,ACTUAL]=counterfactual_intern(obj)

        obj=set_options(obj,varargin{:});
		shocks_db=obj.options.counterfact_shocks_db;
        counterfact_nullify_others=obj.options.counterfact_nullify_others;
        counterfact_regime=obj.options.counterfact_regime;
                
        % get the benchmark history
        [obj,~,~,retcode]=obj.filter();	% by default, this will also smooth....
        
        % extract the state and transition matrices
        if retcode
            error([mfilename,':: model could not be solved'])
        end
        
        shockList={obj.varexo.name};
        endoList=obj.endogenous.name;
        smoothed_shocks=obj.Filters.Expected_smoothed_shocks;
        smoothed_variables=obj.Filters.Expected_smoothed_variables;
        TimeInfo=smoothed_variables.(endoList{1}).TimeInfo;
        if isempty(shocks_db)
            % allow all shocks
            shocks_db=smoothed_shocks;
        end
        if ischar(shocks_db)
            shocks_db=cellstr(shocks_db);
        end
        if isstruct(shocks_db) || isa(shocks_db,'rise_time_series')
            if isa(shocks_db,'rise_time_series')
                shocks_db=pages2struct(shocks_db);
            end
            if ~isequal(shocks_db.(shockList{1}).TimeInfo,smoothed_shocks.(shockList{1}).TimeInfo)
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

        smoothed_probabilities=rise_time_series.collect(obj.Filters.smoothed_probabilities);
        probList=smoothed_probabilities.varnames;
        smoothed_probabilities=permute(double(smoothed_probabilities),[2,1]);
        % apply the smoothed probabilities to SS
        Regimes=obj.Regimes;
        [nreg,nchains]=size(Regimes);
        chainNames={obj.markov_chains.name};
        if isempty(counterfact_regime)
            Probs=ones(nreg,NumberOfObservations);
            for ireg=1:nreg
                reg_i=Regimes(ireg,:);
                for ichain=1:nchains
                    ch_reg=[chainNames{ichain},'_',int2str(reg_i(ichain))];
                    loc=strcmp(ch_reg,probList);
                    Probs(ireg,:)=Probs(ireg,:).*smoothed_probabilities(loc,:);
                end
            end
        else
            if ischar(counterfact_regime)
                counterfact_regime=cellstr(counterfact_regime);
            end
            nn=numel(counterfact_regime)/2;
            if nn~=floor(nn)
                error('counterfact_regime elements must come in pairs')
            end
            counterfact_regime=reshape(counterfact_regime(:)',2,nn);
            choice=true(nreg,1);
            marked=false(1,nchains);
            marked(1)=true; % constant chain always marked
            for kk=1:nn
                ch=counterfact_regime{1,kk};
                if strcmp(ch,'const')
                    continue
                end
                reg=counterfact_regime{2,kk};
                ccol=find(strcmp(ch,chainNames));
                if isempty(ccol)
                    error([ch,' unrecognized as a chain name'])
                end
                if ~ismember(reg,(1:numel(obj.markov_chains(2).states)))
                    error(['regime chosen for chain ',ch,' not in the range of the states of that chain'])
                end
                if marked(ccol)
                    error(['chain ',ch,' cannot be declared more than once'])
                end
                marked(ccol)=true;
                choice=choice & Regimes(:,ccol)==reg;
            end
            not_marked=chainNames(~marked);
            if ~isempty(not_marked)
                disp(not_marked)
                error('the chains listed above have not been declared')
            end
            counterfact_regime=find(choice);
        end
        
        T=obj.T;
        R=obj.R; %
        [endo_nbr,~,solve_expect_order,~]=size(R);
        
        SS=obj.steady_state_and_balanced_growth_path(1:endo_nbr,:);
        % use the new shocks and the smoothed information to reconstruct history
        
        % initialize output
        COUNTER=zeros(endo_nbr,NumberOfObservations);
        % recreate history: things like this could be done in a simulate function
        % but I already have one and this is not the time to mess with it.
        % there is no way to know what y0 is and setting it to be the steady
        % state is only good for filtering, not for smoothing. Hence we impose
        % that the first observation (y1) is given and we compute the rest
        start=2;
        COUNTER(:,1:start-1)=smoothed_variables(:,1:start-1);
        % we then start the iterations at 2
        for t=start:NumberOfObservations
            [A,B,SS_t]=expected_state_matrices(t);
            COUNTER(:,t)=SS_t+A*(COUNTER(:,t-1)-SS_t);
            for r=1:solve_expect_order
                COUNTER(:,t)=COUNTER(:,t)+B(:,:,r)*smoothed_shocks(:,t,r);
            end
        end
        
        COUNTER=double2rise_time_series(COUNTER);
        ACTUAL=double2rise_time_series(smoothed_variables);
        clear smoothed_variables
        
        function [A,B,SS_t]=expected_state_matrices(t)
            if isempty(counterfact_regime)
                probs_t=Probs(:,t);
                A=probs_t(1)*T(:,:,1);
                B=probs_t(1)*R(:,:,:,1);
                SS_t=probs_t(1)*SS(:,1);
                for st=2:nreg
                    A=A+probs_t(st)*T(:,:,st);
                    B=B+probs_t(st)*R(:,:,:,st);
                    SS_t=SS_t+probs_t(st)*SS(:,st);
                end
            else
                A=T(:,:,counterfact_regime);
                B=R(:,:,:,counterfact_regime);
                SS_t=SS(:,counterfact_regime);
            end
        end
        
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



