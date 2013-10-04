function [dsge_irfs,dsge_var_irfs]=irf(obj,varargin)

% 2- If the model is estimated, one may want to draw from the distribution,
% use the mean or the mode for computing the irfs... I should add that too.
% 3- I should add a sub-field for anticipated vs unanticipated irfs?
% 4- shock surprises?

too_small=1e-9;

if isempty(obj)
    dsge_irfs=struct('irf_shock_list','',...
        'irf_var_list','',...
        'irf_periods',40,...
        'irf_anticipate',true,...
        'irf_horizon',1,...
        'irf_param_type','mode',...
        'irf_shock_sign',1,...
        'irf_draws',1,...
        'irf_type','irf',...
        'irf_history',[],...
        'irf_risk',true,...
        'irf_shock_uncertainty',false...
        );%,'irf_parameter_uncertainty',false,'irf_draws',1,...
    return
end
% irf_param_type=obj.options.irf_param_type;

nobj=numel(obj);

dsge_irfs=cell(1,nobj);
dsge_var_irfs=cell(1,nobj);
% check that the models are consistent
check_irf_consistency(obj)
for ii=1:nobj
    [dsge_irfs{ii},dsge_var_irfs{ii}]=irf_intern(obj(ii));
end

dsge_irfs=format_irf_output(dsge_irfs);
dsge_var_irfs=format_irf_output(dsge_var_irfs);

    function [dsge_irfs,dsge_var_irfs]=irf_intern(obj)
        obj=set_options(obj,varargin{:});
        solve_order=obj.options.solve_order;
        
        irf_shock_list        =obj.options.irf_shock_list;
        irf_var_list          =obj.options.irf_var_list  ;
        irf_periods	          =obj.options.irf_periods	  ;
        irf_anticipate        =obj.options.irf_anticipate; % flag to say whether the agents are aware of the shock or not
        irf_horizon           =obj.options.irf_horizon; % this is when the shock hits
        expansion_order       =irf_horizon-1;
        irf_shock_sign        =obj.options.irf_shock_sign;
        irf_type	          =obj.options.irf_type	  ;
        irf_history           =obj.options.irf_history;
        irf_draws	          =obj.options.irf_draws	  ;
        irf_shock_uncertainty	  =obj.options.irf_shock_uncertainty;
        irf_risk	  =obj.options.irf_risk;
        sig=irf_risk;
        
        if isempty(irf_var_list)
            irf_var_list=obj.endogenous.name(obj.endogenous.is_original);% & ~obj.endogenous.is_auxiliary
        elseif ischar(irf_var_list)
            irf_var_list=cellstr(irf_var_list);
        end
        % set the order to the chosen irf_horizon in a way that the
        % shock properties are also affected. But do that only if there is
        % a difference.
        if irf_horizon>obj.options.solve_expect_order
            obj=obj.set_options('solve_expect_order',irf_horizon);
        end
        % otherwise I would need to write 2 lines of code to get the
        % correct results as follows:
        %         obj.options.solve_expect_order=irf_horizon;
        %         [obj.options.shock_properties(:).horizon]=deal(irf_horizon);

        varnames=obj(1).endogenous.name;
        var_locs=locate_variables(irf_var_list,varnames);
        
        h=obj.markov_chains.regimes_number;
        exo_nbr=obj.exogenous.number(1); % unobserved exogenous
        endo_nbr=numel(var_locs);
        % total_endo_nbr=obj.endogenous.number(2);
        
        stochastic=~obj.exogenous.is_observed;
        exoList=obj.exogenous.name(stochastic);
        if isempty(irf_shock_list)
            irf_shock_list=exoList;
        end
        if ischar(irf_shock_list)
            irf_shock_list=cellstr(irf_shock_list);
        end
        position=locate_variables(irf_shock_list,exoList,true);
        if any(isnan(position))
            disp(irf_shock_list(isnan(position)))
            error('The above list of shocks cannot be used in irf')
        end
        nshocks=numel(irf_shock_list);
        
        if ~isempty(irf_history)
            if ~all(ismember(irf_history,1:h))
                error([mfilename,':: all elements in the history of IRFs must be between 1 and ',int2str(h)])
            end
            irf_periods=numel(irf_history);
            irf_type='hirf';
        end
        if h==1
            irf_type='irf';
        end
        girf=strcmp(irf_type,'girf');
        
        % original distribution of states
        %--------------------------------
        % I initialize probabilities in this way in order to
        % allow some variation. The alternative is to start at
        % the ergodic distribution in which case the
        % probabilities will remain unchanged. I would have to
        % use the algorithm hidden in kalman initialization and
        % make it a small function...
        PAI=1/h*ones(h,1);
        
        hstar=1;
        if strcmp(irf_type,'irf'),hstar=h;end
        
        if ismember(irf_type,{'irf','hirf'})
            irf_draws=1;
        end
        if irf_draws==1
            irf_shock_uncertainty=false;
        end
        % initialize
        %-----------
        Impulse_dsge=zeros(endo_nbr,irf_periods+1,nshocks,hstar,irf_draws);
        dsge_var_flag=strcmp(irf_type,'irf') && obj.is_dsge_var_model;
        if dsge_var_flag
            Impulse_dsge_var=zeros(obj.observables.number,irf_periods+1,nshocks,hstar);
        end
        [obj,retcode]=solve(obj);
        if retcode
            error('model cannot be solved')
        end
        % extract the solution and remove the columns corresponding to
        % non-stochastic exogenous variables
        %--------------------------------------------------------------
        solution=obj.solution;
        for state=1:h
            solution.m_e{state}=solution.m_e{state}(:,stochastic,:);
        end
        for gd=1:irf_draws
            for istate=1:hstar
                do_one_dsge_parameter();
                if dsge_var_flag
                    do_one_dsge_var_parameter();
                end
            end
        end
        
        % set to 0 the terms that are too tiny
        Impulse_dsge(abs(Impulse_dsge)<=too_small)=0;
        
        % reshape as time x regimes x variables x shocks x irf_draws
        %------------------------------------------------------------
        % going from variables x time x shocks x regimes x irf_draws
        Impulse_dsge=permute(Impulse_dsge,[2,4,1,3,5]);
        
        % average across  the irf_draws dimension
        %----------------------------------------
        Impulse_dsge=mean(Impulse_dsge,5);
        
        % distribution of irfs
        %---------------------
        startdate=0;
        RegimeNames=strcat('regime_',num2str((1:h)'));
        if ismember(irf_type,{'girf','hirf'})
            RegimeNames=irf_type;
        end
        dsge_irfs=struct();
        dsge_var_irfs=struct();
        for ss=1:exo_nbr
            shock_name=irf_shock_list{ss};
            for vv=1:numel(irf_var_list)
                dsge_irfs.(shock_name).(irf_var_list{vv})=...
                    rise_time_series(startdate,squeeze(Impulse_dsge(:,:,vv,ss)),RegimeNames);
            end
            %----------------------------------------
            if dsge_var_flag
                shock_loc=locate_variables(shock_name,exoList);
                if ss==1
                    % set to 0 the terms that are too tiny
                    Impulse_dsge_var(abs(Impulse_dsge_var)<=too_small)=0;
                    Impulse_dsge_var=permute(Impulse_dsge_var,[2,1,3,4]);% time,varnames,draws,
                end
                for vv_=1:obj.NumberOfObservables(1)
                    dsge_var_irfs.(shock_name).(obj.observables.name{vv_})=...
                        rise_time_series(startdate,squeeze(Impulse_dsge_var(:,vv_,shock_loc,:)));
                end
            end
            %----------------------------------------
        end
        
        function do_one_dsge_var_parameter()
            % impulse responses are computed for all shocks in this
            % case and the parameter draws are outside this function
            param_draw=[];
            Impulse_dsge_var=dsge_var_irf(obj,param_draw,irf_periods,false,false,false);
        end
        function do_one_dsge_parameter()
            for ishock=1:nshocks % user shock list
                if irf_shock_uncertainty
                    e=randn(exo_nbr,expansion_order+1);
                else
                    e=zeros(exo_nbr,expansion_order+1);
                end
                tmp=zeros(nshocks,1);
                this_shock_pos=position(ishock);
                tmp(this_shock_pos)=1*irf_shock_sign;
                e(:,end)=tmp;
                st=choose_state();
                % start at the steady state
                %--------------------------
                x0=solution.ss{st}(:,ones(1,1+girf));
                Impulse_dsge(:,1,this_shock_pos,istate,gd)=x0(var_locs,1);
                for t=1:irf_periods
                    st=choose_state();
                    x0(:,1)=simulation_engine(solution,x0(:,1),e,sig,st,solve_order,irf_anticipate);
                    if girf
                        if t==1
                            e(this_shock_pos,end)=0;
                        end
                        x0(:,2)=simulation_engine(solution,x0(:,2),e,sig,st,solve_order,irf_anticipate);
                        Impulse_dsge(:,t+1,this_shock_pos,istate,gd)=x0(var_locs,1)-x0(var_locs,2);
                        e=randn(exo_nbr,expansion_order+1);
                    else
                        Impulse_dsge(:,t+1,this_shock_pos,istate,gd)=x0(var_locs,1);
                        e=e(:,2:end);
                    end
                end
            end

            function st=choose_state()
                % in the endogenous probability case the configuration of
                % the transition matrix will change ... and this is not
                % implemented here yet.
                if girf
                    % update probabilities
                    %---------------------
                    PAI=solution.Q'*PAI;
                    csp=[0;cumsum(PAI)];
                    st=find(csp>rand,1,'first')-1;
                elseif strcmp(irf_type,'hirf')
                    st=irf_history(t);
                else
                    st=istate;
                end
            end
        end
    end
end

function check_irf_consistency(obj)
nobj=numel(obj);
if nobj>1
    first_list={'endogenous','exogenous'};
    second_list={'regimes'};
    third_list={'irf_shock_list','irf_var_list','irf_periods','irf_type'};
    for ilist=1:numel(first_list)
        ref_list=obj(1).(first_list{ilist}).name;
        for iobj=2:nobj
            list=obj(iobj).(first_list{ilist}).name;
            if ~isequal(ref_list,list)
                warning([first_list{ilist},' is not the same across models'])
            end
        end
    end
    for ilist=1:numel(second_list)
        ref_list=obj(1).markov_chains.(second_list{ilist});
        for iobj=2:nobj
            list=obj(iobj).markov_chains.(second_list{ilist});
            if ~isequal(ref_list,list)
                warning([second_list{ilist},' is not the same across models'])
            end
        end
    end
    for ilist=1:numel(third_list)
        ref_list={obj(1).options.(third_list{ilist})};
        for iobj=2:nobj
            list={obj(iobj).options.(third_list{ilist})};
            if ~isequal(ref_list,list)
                error([third_list{ilist},' should be the same across models'])
            end
        end
    end
end
end

function dsge_irfs=format_irf_output(dsge_irfs)
nobj=numel(dsge_irfs);
if nobj==1
    dsge_irfs=dsge_irfs{1};
else
    if isempty(dsge_irfs{1})
        return
    end
    shockList=fieldnames(dsge_irfs{1});
    tmp=struct();
    for ishock=1:numel(shockList)
        shock_models=cell(1,nobj);
        for mm=1:nobj
            shock_models{mm}=dsge_irfs{mm}.(shockList{ishock});
        end
        tmp.(shockList{ishock})=concatenate_series_from_different_models(shock_models);
    end
    % aggregate
    dsge_irfs=tmp;
end
end

% function [dsge_irfs,dsge_var_irfs]=irf(obj,varargin)
% 
% % 2- If the model is estimated, one may want to draw from the distribution,
% % use the mean or the mode for computing the irfs... I should add that too.
% % 3- I should add a sub-field for anticipated vs unanticipated irfs?
% % 4- shock surprises?
% 
% too_small=1e-9;
% 
% if isempty(obj)
%     dsge_irfs=struct('irf_shock_list','',...
%         'irf_var_list','',...
%         'irf_periods',40,...
%         'irf_anticipate',true,...
%         'irf_horizon',1,...
%         'irf_param_type','mode',...
%         'irf_shock_sign',1,...
%         'irf_draws',1,...
%         'girf_draws',300,...
%         'irf_type','irf',...
%         'irf_history',[],...
%         'irf_risk',true,...
%         'irf_shock_uncertainty',false,...
%         'irf_parameter_uncertainty',false);
%     return
% end
% % irf_param_type=obj.options.irf_param_type;
% 
% nobj=numel(obj);
% 
% dsge_irfs=cell(1,nobj);
% dsge_var_irfs=cell(1,nobj);
% % check that the models are consistent
% check_irf_consistency(obj)
% for ii=1:nobj
%     [dsge_irfs{ii},dsge_var_irfs{ii}]=irf_intern(obj(ii));
% end
% 
% dsge_irfs=format_irf_output(dsge_irfs);
% dsge_var_irfs=format_irf_output(dsge_var_irfs);
% 
%     function [dsge_irfs,dsge_var_irfs]=irf_intern(obj)
%         obj=set_options(obj,varargin{:});
%         solve_order=obj.options.solve_order;
%         
%         irf_shock_list        =obj.options.irf_shock_list;
%         irf_var_list          =obj.options.irf_var_list  ;
%         irf_periods	          =obj.options.irf_periods	  ;
%         irf_anticipate        =obj.options.irf_anticipate; % flag to say whether the agents are aware of the shock or not
%         irf_horizon           =obj.options.irf_horizon; % this is when the shock hits
%         expansion_order       =irf_horizon-1;
%         irf_shock_sign        =obj.options.irf_shock_sign;
%         % how many parameter draws
%         irf_draws	          =obj.options.irf_draws	  ;
%         irf_type	          =obj.options.irf_type	  ;
%         irf_history           =obj.options.irf_history;
%         girf_draws	          =obj.options.girf_draws	  ;
%         irf_shock_uncertainty	  =obj.options.irf_shock_uncertainty;
%         irf_parameter_uncertainty =obj.options.irf_parameter_uncertainty;
%         irf_risk	  =obj.options.irf_risk;
%         sig=irf_risk;
%         
%         if isempty(irf_var_list)
%             irf_var_list=obj.endogenous.name(obj.endogenous.is_original);% & ~obj.endogenous.is_auxiliary
%         elseif ischar(irf_var_list)
%             irf_var_list=cellstr(irf_var_list);
%         end
%         % set the order to the chosen irf_horizon in a way that the
%         % shock properties are also affected. But do that only if there is
%         % a difference.
%         if irf_horizon>obj.options.solve_expect_order
%             obj=obj.set_options('solve_expect_order',irf_horizon);
%         end
%         % otherwise I would need to write 2 lines of code to get the
%         % correct results as follows:
%         %         obj.options.solve_expect_order=irf_horizon;
%         %         [obj.options.shock_properties(:).horizon]=deal(irf_horizon);
% 
%         varnames=obj(1).endogenous.name;
%         var_locs=locate_variables(irf_var_list,varnames);
%         
%         h=obj.markov_chains.regimes_number;
%         exo_nbr=obj.exogenous.number(1); % unobserved exogenous
%         endo_nbr=numel(var_locs);
%         % total_endo_nbr=obj.endogenous.number(2);
%         
%         stochastic=~obj.exogenous.is_observed;
%         exoList=obj.exogenous.name(stochastic);
%         if isempty(irf_shock_list)
%             irf_shock_list=exoList;
%         end
%         if ischar(irf_shock_list)
%             irf_shock_list=cellstr(irf_shock_list);
%         end
%         position=locate_variables(irf_shock_list,exoList,true);
%         if any(isnan(position))
%             disp(irf_shock_list(isnan(position)))
%             error('The above list of shocks cannot be used in irf')
%         end
%         nshocks=numel(irf_shock_list);
%         
%         if ~isempty(irf_history)
%             if ~all(ismember(irf_history,1:h))
%                 error([mfilename,':: all elements in the history of IRFs must be between 1 and ',int2str(h)])
%             end
%             irf_periods=numel(irf_history);
%             irf_type='hirf';
%         end
%         if h==1
%             irf_type='irf';
%         end
%         girf=strcmp(irf_type,'girf');
%         
%         % original distribution of states
%         %--------------------------------
%         % I initialize probabilities in this way in order to
%         % allow some variation. The alternative is to start at
%         % the ergodic distribution in which case the
%         % probabilities will remain unchanged. I would have to
%         % use the algorithm hidden in kalman initialization and
%         % make it a small function...
%         PAI=1/h*ones(h,1);
%         
%         hstar=1;
%         if strcmp(irf_type,'irf'),hstar=h;end
%         
%         if ismember(irf_type,{'irf','hirf'})
%             girf_draws=1;
%         end
%         
%         % initialize
%         %-----------
%         Impulse_dsge=zeros(endo_nbr,irf_periods+1,nshocks,hstar,irf_draws,girf_draws);
%         dsge_var_flag=strcmp(irf_type,'irf') && obj.is_dsge_var_model;
%         if dsge_var_flag
%             Impulse_dsge_var=zeros(obj.observables.number,irf_periods+1,nshocks,hstar,irf_draws);
%         end
%         d=0;
%         while d<irf_draws
%             d=d+1;
%             % draw a parameter and solve the model
%             %-------------------------------------
%             if irf_parameter_uncertainty
%                 param_draw=draw_parameter();
%                 % push the parameter vector into the dsge model
%                 %----------------------------------------------
%                 obj=assign_estimates(obj,param_draw);
%             end
%             [obj,retcode]=solve(obj);
%             if retcode
%                 if ~irf_parameter_uncertainty
%                     error('model cannot be solved')
%                 end
%                 d=d-1;
%                 continue
%             end
%             % extract the solution and remove the columns corresponding to
%             % non-stochastic exogenous variables
%             %--------------------------------------------------------------
%             solution=obj.solution;
%             for state=1:h
%                 solution.m_e{state}=solution.m_e{state}(:,stochastic,:);
%             end
%             for gd=1:girf_draws
%                 for istate=1:hstar
%                     do_one_dsge_parameter();
%                     if dsge_var_flag
%                         do_one_dsge_var_parameter();
%                     end
%                 end
%             end
%         end
%         
%         % set to 0 the terms that are too tiny
%         Impulse_dsge(abs(Impulse_dsge)<=too_small)=0;
%         
%         % reshape as time x regimes x variables x shocks x irf_draws x girf_draws
%         %------------------------------------------------------------------------
%         % going from variables x time x shocks x regimes x irf_draws x girf_draws
%         Impulse_dsge=permute(Impulse_dsge,[2,4,1,3,5,6]);
%         
%         % remove the girf dimension
%         %--------------------------
%         Impulse_dsge=mean(Impulse_dsge,6);
%         
%         % sort along the parameter drawing dimension
%         %-------------------------------------------
%         Impulse_dsge=sort(Impulse_dsge,5);
%         
%         % % average across all draws for the mean
%         % %--------------------------------------
%         % Imp_mean=mean(Imp,5);
%         
%         % distribution of irfs
%         %---------------------
%         startdate=0;
%         RegimeNames=strcat('regime_',num2str((1:h)'));
%         if ismember(irf_type,{'girf','hirf'})
%             RegimeNames=irf_type;
%         end
%         dsge_irfs=struct();
%         dsge_var_irfs=struct();
%         for ss=1:exo_nbr
%             shock_name=irf_shock_list{ss};
%             for vv=1:numel(irf_var_list)
%                 dsge_irfs.(shock_name).(irf_var_list{vv})=...
%                     rise_time_series(startdate,squeeze(Impulse_dsge(:,:,vv,ss,:)),RegimeNames);
%             end
%             %----------------------------------------
%             if dsge_var_flag
%                 shock_loc=locate_variables(shock_name,exoList);
%                 if ss==1
%                     % set to 0 the terms that are too tiny
%                     Impulse_dsge_var(abs(Impulse_dsge_var)<=too_small)=0;
%                     Impulse_dsge_var=permute(Impulse_dsge_var,[2,1,3,4]);% time,varnames,draws,
%                 end
%                 for vv_=1:obj.NumberOfObservables(1)
%                     dsge_var_irfs.(shock_name).(obj.observables.name{vv_})=...
%                         rise_time_series(startdate,squeeze(Impulse_dsge_var(:,vv_,:,shock_loc)));
%                 end
%             end
%             %----------------------------------------
%         end
%         
%         function do_one_dsge_var_parameter()
%             % impulse responses are computed for all shocks in this
%             % case
%             Impulse_dsge_var(:,:,dd,:)=dsge_var_irf(obj,param_draw,irf_periods,false,false,false);
%         end
%         function do_one_dsge_parameter()
%             for ishock=1:nshocks % user shock list
%                 if irf_shock_uncertainty
%                     e=randn(exo_nbr,expansion_order+1);
%                 else
%                     e=zeros(exo_nbr,expansion_order+1);
%                 end
%                 tmp=zeros(nshocks,1);
%                 this_shock_pos=position(ishock);
%                 tmp(this_shock_pos)=1*irf_shock_sign;
%                 e(:,end)=tmp;
%                 st=choose_state();
%                 % start at the steady state
%                 %--------------------------
%                 x0=solution.ss{st}(:,ones(1,1+girf));
%                 Impulse_dsge(:,1,this_shock_pos,istate,d,gd)=x0(var_locs,1);
%                 for t=1:irf_periods
%                     st=choose_state();
%                     x0(:,1)=simulation_engine(solution,x0(:,1),e,sig,st,solve_order,irf_anticipate);
%                     if girf
%                         if t==1
%                             e(this_shock_pos,end)=0;
%                         end
%                         x0(:,2)=simulation_engine(solution,x0(:,2),e,sig,st,solve_order,irf_anticipate);
%                         Impulse_dsge(:,t+1,this_shock_pos,istate,d,gd)=x0(var_locs,1)-x0(var_locs,2);
%                         e=randn(exo_nbr,expansion_order+1);
%                     else
%                         Impulse_dsge(:,t+1,this_shock_pos,istate,d,gd)=x0(var_locs,1);
%                         e=e(:,2:end);
%                     end
%                 end
%             end
% 
%             function st=choose_state()
%                 % in the endogenous probability case the configuration of
%                 % the transition matrix will change ... and this is not
%                 % implemented here yet.
%                 if girf
%                     % update probabilities
%                     %---------------------
%                     PAI=solution.Q'*PAI;
%                     csp=[0;cumsum(PAI)];
%                     st=find(csp>rand,1,'first')-1;
%                 elseif strcmp(irf_type,'hirf')
%                     st=irf_history(t);
%                 else
%                     st=istate;
%                 end
%             end
%         end
%     end
% end
% 
% function check_irf_consistency(obj)
% nobj=numel(obj);
% if nobj>1
%     first_list={'endogenous','exogenous'};
%     second_list={'regimes'};
%     third_list={'irf_shock_list','irf_var_list','irf_periods','irf_type'};
%     for ilist=1:numel(first_list)
%         ref_list=obj(1).(first_list{ilist}).name;
%         for iobj=2:nobj
%             list=obj(iobj).(first_list{ilist}).name;
%             if ~isequal(ref_list,list)
%                 warning([first_list{ilist},' is not the same across models'])
%             end
%         end
%     end
%     for ilist=1:numel(second_list)
%         ref_list=obj(1).markov_chains.(second_list{ilist});
%         for iobj=2:nobj
%             list=obj(iobj).markov_chains.(second_list{ilist});
%             if ~isequal(ref_list,list)
%                 warning([second_list{ilist},' is not the same across models'])
%             end
%         end
%     end
%     for ilist=1:numel(third_list)
%         ref_list={obj(1).options.(third_list{ilist})};
%         for iobj=2:nobj
%             list={obj(iobj).options.(third_list{ilist})};
%             if ~isequal(ref_list,list)
%                 error([third_list{ilist},' should be the same across models'])
%             end
%         end
%     end
% end
% end
% 
% function dsge_irfs=format_irf_output(dsge_irfs)
% nobj=numel(dsge_irfs);
% if nobj==1
%     dsge_irfs=dsge_irfs{1};
% else
%     if isempty(dsge_irfs{1})
%         return
%     end
%     shockList=fieldnames(dsge_irfs{1});
%     tmp=struct();
%     for ishock=1:numel(shockList)
%         shock_models=cell(1,nobj);
%         for mm=1:nobj
%             shock_models{mm}=dsge_irfs{mm}.(shockList{ishock});
%         end
%         tmp.(shockList{ishock})=concatenate_series_from_different_models(shock_models);
%     end
%     % aggregate
%     dsge_irfs=tmp;
% end
% end