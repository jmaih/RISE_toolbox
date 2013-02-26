function [dsge_irfs,dsge_var_irfs]=irf(obj,varargin)

% 2- If the model is estimated, one may want to draw from the distribution,
% use the mean or the mode for computing the irfs... I should add that too.
% 3- I should add a sub-field for anticipated vs unanticipated irfs?
% 4- shock surprises?

too_small=1e-9;

Defaults=struct('irf_shock_list','',...
    'irf_var_list','',...
    'irf_periods',40,...
    'irf_anticipate',true,...
    'irf_horizon',1,...
    'irf_param_type','mode',...
    'irf_shock_sign',1,...
    'irf_draws',1,...
    'girf_draws',300,...
    'irf_type','irf',...
    'irf_history',[]);

if isempty(obj)
    dsge_irfs=Defaults;
    return
end
nobj=numel(obj);
Fields=fieldnames(Defaults);

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
%         obj.options=mysetfield(obj.options,varargin{:});
        oldoptions=obj.options;
        thisDefault=Defaults;
        for ifield=1:numel(Fields)
            vi=Fields{ifield};
            thisDefault.(vi)=oldoptions.(vi);
        end
        irf_shock_list=thisDefault.irf_shock_list;
        irf_var_list  =thisDefault.irf_var_list  ;
        irf_periods	  =thisDefault.irf_periods	  ;
        irf_anticipate=thisDefault.irf_anticipate;
        irf_horizon=thisDefault.irf_horizon;
        irf_param_type=thisDefault.irf_param_type;
        irf_shock_sign=thisDefault.irf_shock_sign;
        % how many parameter draws
        irf_draws	  =thisDefault.irf_draws	  ;
        girf_draws	  =thisDefault.girf_draws	  ;
        irf_type	  =thisDefault.irf_type	  ;
        irf_history   =thisDefault.irf_history;
        if isempty(irf_shock_list)
            irf_shock_list={obj(1).varexo.name};
        elseif ischar(irf_shock_list)
            irf_shock_list=cellstr(irf_shock_list);
        end
        if isempty(irf_var_list)
            irf_var_list={obj.orig_varendo.name};
        elseif ischar(irf_var_list)
            irf_var_list=cellstr(irf_var_list);
        end
        shocks_id=locate_variables(irf_shock_list,{obj(1).varexo.name});
        
        varnames={obj(1).varendo.name};
        var_locs=locate_variables(irf_var_list,varnames);
        
        shock_only_one=false;
%         relative=false;
        
        % set the order to the chosen irf_horizon in a way that the
        % shock properties are also affected. But do that only if there is
        % a difference.
        if irf_horizon>obj.options.order
            obj=obj.set_options('order',irf_horizon);
        end
        % otherwise I would need to write 2 lines of code to get the
        % correct results as follows:
        %         obj.options.order=irf_horizon;
        %         [obj.options.shock_properties(:).horizon]=deal(irf_horizon);
        
        % number of endogenous variables after solving...
        NumberOfEndogenous=obj.NumberOfEndogenous(2);
        NumberOfShocks=numel(shocks_id);
        NumberOfRegimes=obj.NumberOfRegimes;
        if ~isempty(irf_history)
            if ~all(ismember(irf_history,1:NumberOfRegimes))
                error([mfilename,':: all elements in the history of IRFs must be between 1 and ',int2str(NumberOfRegimes)])
            end
            irf_periods=numel(irf_history);
            irf_type='hirf';
        end
        
        if NumberOfRegimes==1
            irf_type='irf';
        end
        d=0;
        dsge_var_irfs=rise_time_series.empty(0,1);
        switch irf_type
            case 'irf'
                Z=zeros(NumberOfEndogenous,irf_periods+1,irf_draws,NumberOfShocks,NumberOfRegimes);
                if obj.is_dsge_var_model
                    Z_bvar_dsge=zeros(obj.NumberOfObservables,irf_periods+1,irf_draws,NumberOfShocks);
                end
            case {'girf','hirf'}
                Z=zeros(NumberOfEndogenous,irf_periods+1,irf_draws,NumberOfShocks);
            otherwise
                error([mfilename,':: unknown option for impulse response(',obj.options.irf_type,...
                    '). Valid options are ''irf'', ''girf'' and ''hirf'''])
        end
        while d<irf_draws
            if d==0
                obj=obj.set_parameters('mode');
                param_draw=[];
            else
                % this function has not been written yet. But how do I do the
                % metropolis draws, etc. do I save them?
                [obj,param_draw]=obj.draw_parameters(irf_param_type);
            end
            [obj,retcode]=solve(obj);
            if retcode
                if irf_draws==1
                    error([mfilename,':: model could not be solved'])
                else
                    continue
                end
            end
            d=d+1;
            
%             steady_state=vertcat(obj.varendo.det_steady_state);
%             same_steady_state=max(max(abs(bsxfun(@minus,steady_state,steady_state(:,1)))))<1e-9;
            
            T=obj.T;
            R=obj.R;
            Q=obj.Q;
            for ss=1:NumberOfShocks
                %                 shock_name=irf_shock_list{ss}; %#ok<USENS>
                shock_loc=shocks_id(ss);
                switch irf_type
                    case 'girf'
                        for dd=1:girf_draws
                            y0=[];
                            PAI=1/NumberOfRegimes*ones(NumberOfRegimes,1);
                            % I initialize probabilities in this way in order to
                            % allow some variation. The alternative is to start at
                            % the ergodic distribution in which case the
                            % probabilities will remain unchanged. I would have to
                            % use the algorithm hidden in kalman initialization and
                            % make it a small function...
                            for t=1:irf_periods
                                [imp,y0]=update_impulse(y0,T,R,PAI,shock_loc,t,irf_shock_sign,shock_only_one);
                                Z(:,t+1,d,ss)=Z(:,t+1,d,ss)+imp;
                                % update probabilities
                                PAI=Q'*PAI;
                            end
                        end
                        Z(:,:,d,ss)=Z(:,:,d,ss)/girf_draws;
% % %                         if same_steady_state
% % %                             Z(:,:,d,ss)=bsxfun(@plus,Z(:,:,d,ss),relative*steady_state(:,1));
% % %                         else
% % %                             disp([mfilename,':: dealing with switching steady states is an unsolved problem here for the moment'])
% % %                             % I don't know how to deal with multiple steady states
% % %                             % in this case
% % %                         end
                    case 'irf'
                        for reg=1:NumberOfRegimes
                            % create one page for each of the shocks in each regime
                            y0=[];
                            PAI=1;
                            for t=1:irf_periods
                                [Z(:,t+1,d,ss,reg),y0]=update_impulse(y0,T(:,:,reg),R(:,:,:,reg),PAI,...
                                    shock_loc,t,irf_shock_sign);
                            end
% % %                             Z(:,:,d,ss,reg)=bsxfun(@plus,Z(:,:,d,ss,reg),relative*steady_state(:,reg));
                        end
                    case 'hirf'
                        % create one page for each of the shocks in each regime
                        y0=[];
                        PAI=1;
                        for t=1:irf_periods
                            reg=irf_history(t);
                            [Z(:,t+1,d,ss),y0]=update_impulse(y0,T(:,:,reg),R(:,:,:,reg),PAI,...
                                shock_loc,t,irf_shock_sign);
                        end
                end
            end
            if strcmp(irf_type,'irf') && obj.is_dsge_var_model
                % impulse responses are computed for all shocks in this
                % case
                Z_bvar_dsge(:,:,dd,:)=dsge_var_irf(obj,param_draw,irf_periods,false,false,false);
                %                 Z_bvar_dsge=zeros(obj.NumberOfObservables,irf_periods+1,irf_draws,NumberOfShocks);
            end
        end
        
        Z=permute(Z,[2,5,3,1,4]); % time,regimes,draws,varnames,shocks
        % set to 0 the terms that are too tiny
        Z(abs(Z)<=too_small)=0;
        startdate=0;
        RegimeNames=strcat('regime_',num2str((1:NumberOfRegimes)'));
        if ismember(irf_type,{'girf','hirf'})
            RegimeNames=irf_type;
        end
        for ss=1:NumberOfShocks
            shock_name=irf_shock_list{ss};
            for vv=1:numel(irf_var_list)
                dsge_irfs.(shock_name).(irf_var_list{vv})=...
                    rise_time_series(startdate,Z(:,:,:,var_locs(vv),ss),RegimeNames);
                if ss==1 && vv==1
                    startdate=dsge_irfs.(shock_name).(irf_var_list{vv}).TimeInfo;
                end
            end
            shock_loc=shocks_id(ss);
            if strcmp(irf_type,'irf') && obj.is_dsge_var_model
                if ss==1
                    Z_bvar_dsge=permute(Z_bvar_dsge,[2,1,3,4]);% time,varnames,draws,
                end
                disp([mfilename,':: check the syntax here for the number observables and the names'])
                % set to 0 the terms that are too tiny
                Z_bvar_dsge(abs(Z_bvar_dsge)<=too_small)=0;
                for vv=1:obj.NumberOfObservables(1)
                    dsge_var_irfs.(shock_name).(obj.varobs(vv).name)=...
                        rise_time_series(startdate,squeeze(Z_bvar_dsge(:,vv,:,shock_loc)));
                end
            end
        end
        clear Z Z_bvar_dsge
        
        function [imp,y0]=update_impulse(y0,T,R,PAI,shock_id,t,irf_shock_sign,shock_only_one)
            if nargin<8
                shock_only_one=true;
                if nargin<7
                    irf_shock_sign=1;
                end
            end
            [n,exo_nbr,~,h]=size(R);
            if h>1
                % pick a state
                csp=[0;cumsum(PAI)];
                state=find(csp>rand,1,'first')-1;
            else
                state=1;
            end
            if isempty(y0)
                y0=zeros(n,1);
                if h>1
                    y0=[y0,y0];
                end
                % 1 standard-deviation shock. That is unanticipated
                % the first time it is known but which will remain
                % in the system
                if irf_anticipate || irf_horizon==1
                    y0(:,1)=irf_shock_sign*R(:,shock_id,irf_horizon,state);
                end
            else
                y0(:,1)=T(:,:,state)*y0(:,1);
                if h>1 % then we are doing GIRF
                    y0(:,2)=T(:,:,state)*y0(:,2);
                    % agents get hit in the current period by new
                    % shocks that are unanticipated. I don't know
                    % whether I should include all other possible
                    % shocks...
                    if shock_only_one
                        shock=randn;
                        y0(:,1)=y0(:,1)+R(:,shock_id,1,state)*shock;
                        y0(:,2)=y0(:,2)+R(:,shock_id,1,state)*shock;
                    else
                        shock=randn(exo_nbr,1);
                        y0(:,1)=y0(:,1)+R(:,:,1,state)*shock;
                        y0(:,2)=y0(:,2)+R(:,:,1,state)*shock;
                    end
                end
                if t<=irf_horizon
                    step=irf_horizon-t+1;
                    if irf_anticipate || step==1
                        y0(:,1)=y0(:,1)+irf_shock_sign*R(:,shock_id,step,state);
                    end
                end
            end
            imp=y0(:,1);
            if h>1
                imp=imp-y0(:,2);
            end
        end
    end

end

function check_irf_consistency(obj)
nobj=numel(obj);
if nobj>1
    first_list={'varendo','varexo','Regimes','markov_chains'};
    second_list={'irf_shock_list','irf_var_list','irf_periods','irf_type'};
    for ilist=1:numel(first_list)
        ref_list={obj(1).(first_list{ilist})};
        for iobj=2:nobj
            list={obj(iobj).(first_list{ilist})};
            if ~isequal(ref_list,list)
                warning([first_list{ilist},' is not the same across models'])
            end
        end
    end
    for ilist=1:numel(second_list)
        ref_list={obj(1).options.(second_list{ilist})};
        for iobj=2:nobj
            list={obj(iobj).options.(second_list{ilist})};
            if ~isequal(ref_list,list)
                error([second_list{ilist},' should be the same across models'])
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

% function dsge_irfs=format_irf_output(dsge_irfs)
% nobj=numel(dsge_irfs);
% if nobj==1
%     dsge_irfs=dsge_irfs{1};
% else
%     if isempty(dsge_irfs{1})
%         return
%     end
%     tmp=struct();
%     shockList=fieldnames(dsge_irfs{1});
%     varList=fieldnames(dsge_irfs{1}.(shockList{1}));
%     mydates=dsge_irfs{1}.(shockList{1}).(varList{1}).TimeInfo;
%     for ishock=1:numel(shockList)
%         for ivar=1:numel(varList)
%             for imod=1:nobj
%                 datta=double(dsge_irfs{imod}.(shockList{ishock}).(varList{ivar}));
%                 if imod==1 && ivar==1 && ishock==1
%                     irf_size=size(datta);
%                     tank=zeros(irf_size(1),nobj,irf_size(2));
%                     mod_names=strcat('model_',cellfun(@num2str,num2cell(1:nobj),'uniformOutput',false));
%                     reg_names=strcat('regime_',cellfun(@num2str,num2cell(1:irf_size(2)),'uniformOutput',false));
%                 end
%                 tank(:,imod,:)=datta;
%             end
%             if irf_size(2)>1
%                 for ireg=1:irf_size(2)
%                     tmp.(shockList{ishock}).(reg_names{ireg}).(varList{ivar})=...
%                         rise_time_series(mydates,tank(:,:,ireg),mod_names);
%                 end
%             else
%                 tmp.(shockList{ishock}).(varList{ivar})=...
%                     rise_time_series(mydates,tank,mod_names);
%             end
%         end
%     end
%     % aggregate
%     dsge_irfs=tmp;
% end
% 
% end
