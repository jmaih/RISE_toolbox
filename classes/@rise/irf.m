function [dsge_irfs,dsge_var_irfs]=irf(obj,varargin)

% 2- If the model is estimated, one may want to draw from the distribution,
% use the mean or the mode for computing the irfs... I should add that too.
% 3- I should add a sub-field for anticipated vs unanticipated irfs?
% 4- shock surprises?

Defaults=struct('irf_shock_list','',...
    'irf_var_list','',...
    'irf_periods',40,...
    'irf_anticipate',1,...
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
for ii=1:nobj
    [dsge_irfs{ii},dsge_var_irfs{ii}]=irf_intern(obj(ii));
end


    function [dsge_irfs,dsge_var_irfs]=irf_intern(obj)
        
        oldoptions=obj.options;
        oldoptions=mysetfield(oldoptions,varargin{:});
        thisDefault=Defaults;
        for ifield=1:numel(Fields)
            vi=Fields{ifield};
            thisDefault.(vi)=oldoptions.(vi);
        end
        irf_shock_list=thisDefault.irf_shock_list;
        irf_var_list  =thisDefault.irf_var_list  ;
        irf_periods	  =thisDefault.irf_periods	  ;
        irf_anticipate=thisDefault.irf_anticipate;
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
        relative=false;
        
        % set the order to the chosen irf_anticipate in a way that the
        % shock properties are also affected. But do that only if there is
        % a difference.
        if irf_anticipate>obj.options.order
            obj=obj.set_options('order',irf_anticipate);
        end
        % otherwise I would need to write 2 lines of code to get the
        % correct results as follows:
        %         obj.options.order=irf_anticipate;
        %         [obj.options.shock_properties(:).horizon]=deal(irf_anticipate);
        
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
            
            steady_state=vertcat(obj.varendo.det_steady_state);
            same_steady_state=max(max(abs(bsxfun(@minus,steady_state,steady_state(:,1)))))<1e-9;
            
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
                                [imp,y0]=update_impulse(y0,T,R,PAI,irf_anticipate,shock_loc,t,irf_shock_sign,shock_only_one);
                                Z(:,t+1,d,ss)=Z(:,t+1,d,ss)+imp;
                                % update probabilities
                                PAI=Q'*PAI;
                            end
                        end
                        Z(:,:,d,ss)=Z(:,:,d,ss)/girf_draws;
                        if same_steady_state
                            Z(:,:,d,ss)=bsxfun(@plus,Z(:,:,d,ss),relative*steady_state(:,1));
                        else
                            disp([mfilename,':: dealing with switching steady states is an unsolved problem here for the moment'])
                            % I don't know how to deal with multiple steady states
                            % in this case
                        end
                    case 'irf'
                        for reg=1:NumberOfRegimes
                            % create one page for each of the shocks in each regime
                            y0=[];
                            PAI=1;
                            for t=1:irf_periods
                                [Z(:,t+1,d,ss,reg),y0]=update_impulse(y0,T(:,:,reg),R(:,:,:,reg),PAI,...
                                    irf_anticipate,shock_loc,t,irf_shock_sign);
                            end
                            Z(:,:,d,ss,reg)=bsxfun(@plus,Z(:,:,d,ss,reg),relative*steady_state(:,reg));
                        end
                    case 'hirf'
                        % create one page for each of the shocks in each regime
                        y0=[];
                        PAI=1;
                        for t=1:irf_periods
                            reg=irf_history(t);
                            [Z(:,t+1,d,ss),y0]=update_impulse(y0,T(:,:,reg),R(:,:,:,reg),PAI,...
                                irf_anticipate,shock_loc,t,irf_shock_sign);
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
                for vv=1:obj.NumberOfObservables(1)
                    dsge_var_irfs.(shock_name).(obj.varobs(vv).name)=...
                        rise_time_series(startdate,squeeze(Z_bvar_dsge(:,vv,:,shock_loc)));
                end
            end
        end
        clear Z Z_bvar_dsge
    end
end

function [imp,y0]=update_impulse(y0,T,R,PAI,irf_anticipate,shock_id,t,irf_shock_sign,shock_only_one)
if nargin<9
    shock_only_one=true;
    if nargin<8
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
    y0(:,1)=irf_shock_sign*R(:,shock_id,irf_anticipate,state)*1;
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
    if t<=irf_anticipate
        y0(:,1)=y0(:,1)+irf_shock_sign*R(:,shock_id,irf_anticipate-t+1,state);
        if h>1
            y0(:,2)=y0(:,2)+irf_shock_sign*R(:,shock_id,irf_anticipate-t+1,state);
        end
    end
end
imp=y0(:,1);
if h>1
    imp=imp-y0(:,2);
end
end

% r0=obj(1).options.graphics(1);
% c0=obj(1).options.graphics(2);
% selected_var_tex_list={obj(1).varendo(var_locs).tex_name};
% nfields=nobj;
% for jj=1:numel(shocks_id)
%     shock_name=irf_shock_list{jj};
%     X=[];
%     for ii=1:nobj
%         if ii==1 && jj==1 % check the kind of irfs
%             if isfield(obj(1).irfs,'regime_girf')
%                 string='irfs.regime_girf';
%                 regime=1;
%             elseif isfield(obj(1).irfs,'regime_1')
%                 string='irfs.regime_';
%                 nfields=numel(fieldnames(obj(1).irfs));
%                 if nfields>1 && nobj>1
%                     return
%                 end
%                 regime=num2str(transpose(1:nfields));
%             elseif isfield(obj(1).dsge_var,'irfs')
%                 string='dsge_var.irfs';
%             end
%         end
%         if strcmp(string,'irfs.regime_')
%             tmp=[];
%             for kk=1:nfields
%                 eval(['tmp1=double(obj(ii).',string,int2str(kk),'.(shock_name));'])
%                 tmp=[tmp,mean(tmp1(:,var_locs,:),3)]; %#ok<NODEF,*AGROW>
%             end
%         else
%             eval(['tmp=double(obj(ii).',string,'.(shock_name));'])
%             tmp=mean(tmp(:,var_locs,:),3);
%         end
%         X=[X,tmp];
%     end
%     X=reshape(X,[irf_periods+1,numel(var_locs),nfields]);
%     %     plot_impulse_response(X,selected_var_tex_list,shock_name,regime,r0,c0,irf_anticipate)
% end
% rows and columns in graphs probably not needed?


% function plot_impulse_response(X,var_tex_list,shock_name,regime,r0,c0,irf_anticipate)
%
% nn=numel(var_tex_list);
% nstar=r0*c0;
% nfig=ceil(nn/nstar);
% X(abs(X)<1e-9)=0;
%
% titel=['Orthogonalized shocks to ',shock_name];
% for fig=1:nfig
%     if nfig>1
%         titelfig=[titel,' ',int2str(fig)];
%     else
%         titelfig=titel;
%     end
%     figure('name',titelfig);
%     [Remains,r,c]=number_of_rows_and_columns_in_figure(fig,nn,r0,c0);
%     for ii=1:min(nstar,Remains)
%         var_id=(fig-1)*nstar+ii;
%         data=squeeze(X(:,var_id,:));
%         subplot(r,c,ii)
%         hold on
%         plot(data,'linewidth',2)
%         title(var_tex_list{var_id},'fontsize',12)%,'interpreter','none'
%         if ii==1 && numel(regime)>1
%             legend(regime)
%         end
%         plot([1,length(data)],[data(1),data(1)],'r','linewidth',2) % plot([0,irf_anticipate],[min(min(data)),max(max(data))],'r')
%         %         plot([0,irf_anticipate],[0,0],'r','linewidth',2) % plot([0,irf_anticipate],[min(min(data)),max(max(data))],'r')
%         hold off
%         axis tight
%     end
% end
% end



