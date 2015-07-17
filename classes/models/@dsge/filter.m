function [obj,LogLik,Incr,retcode]=filter(obj,varargin)
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
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    obj=filter_initialization(obj);
    return
end

nobj=numel(obj);
Incr=[];
if nobj>1
    retcode=nan(1,nobj);
    LogLik=nan(1,nobj);
    for iobj=1:nobj
        [obj(iobj),LogLik(iobj),~,retcode(iobj)]=filter(obj(iobj),varargin{:});
    end
    return
end
% initialize remaining outputs
%-------------------------------
LogLik=-obj.options.estim_penalty;

if ~isempty(varargin)
    obj=set(obj,varargin{:});
end

if ~obj.data_are_loaded
    obj=obj.load_data;
end

%% solve the object
if obj.warrant_resolving
    [obj,retcode]=obj.solve();
    if retcode
        if obj(1).options.debug
            utils.error.decipher(retcode)
        end
        return
    end
end
h=obj.markov_chains.regimes_number;

%% initialize: load solution, set initial conditions, constraints, etc.
[init,retcode]=filter_initialization(obj);

if retcode
    if obj(1).options.debug
        utils.error.decipher(retcode)
    end
    return
end
%% re-inflate everything in order to store filters if necessary
init=re_inflator(init,obj.options.kf_filtering_level);

%% extract data and update the position of observables
%------------------------------------------------------
data_info=obj.data;
data_info.varobs_id=obj.inv_order_var(data_info.varobs_id);

data_structure=data_info.data_structure;
% no_more_missing=data_info.no_more_missing;
% over-write data_info.varobs_id
data_info.varobs_id=init.obs_id;
obs_id=data_info.varobs_id;
y=data_info.y;
% exogenous data will have only one page...
%--------------------------------------------
U=data_info.x(:,:,1);

% state matrices
%---------------
init.ff=do_one_step_forecast(init.T,init.steady_state,init.sep_compl,...
    init.anticipated_shocks,init.state_vars_location,obj.options.simul_sig,...
    init.horizon,init.is_det_shock);%,shoot

% mapping from states to observables
%------------------------------------
z=recover_positions(obs_id,data_structure,data_info.no_more_missing,...
    data_info.last_good_conditional_observation);
is_real_time=isfield(data_info,'npages') && data_info.npages>1 &&...
    any([~isempty(data_info.restr_y_id),~isempty(data_info.restr_z_id)]);
if isempty(obj.options.kf_user_algo)
    if obj.options.solve_order==1
        if is_real_time
            if ~isempty(U)
                error('real-time filtering with exogenous variables not ready')
            end
            mu_id=data_info.restr_y_id;
            iov=obj.inv_order_var;
            mu_id=iov(real(mu_id))+imag(mu_id)*1i;
            shock_id=data_info.restr_z_id;
            e_data={data_info.z,shock_id};
            y_data={y,mu_id};
            [LogLik,Incr,retcode,Filters]=msre_kalman_cell_real_time(...
                init,y_data,U,z,e_data,obj.options);
        else
            [LogLik,Incr,retcode,Filters]=constrained_regime_switching_kalman_filter_cell(...
                init,y,U,z,obj.options);
        end
    else
        if is_real_time
            error('nonlinear filtering with real-time information not ready')
        end
        [LogLik,Incr,retcode,Filters]=switching_divided_difference_filter(...
            init,y,U,z,obj.options);
    end
else
    vargs={};
    if ischar(obj.options.kf_user_algo)
        obj.options.kf_user_algo=str2func(obj.options.kf_user_algo);
    end
    if iscell(obj.options.kf_user_algo)
        if ischar(obj.options.kf_user_algo{1})
            obj.options.kf_user_algo{1}=str2func(obj.options.kf_user_algo{1});
        end
        vargs=obj.options.kf_user_algo(2:end);
        user_filter=obj.options.kf_user_algo{1};
    else
        user_filter=obj.options.kf_user_algo;
    end
    [LogLik,Incr,retcode,Filters]=user_filter(init,y,U,z,obj.options,vargs{:});
end

if obj.options.kf_filtering_level && ~retcode
    % squash the filters from the inflation
    %---------------------------------------
    table_map={
        'a','P','eta_tlag'
        'att','Ptt','eta_tt'
        'atT','PtT','eta'
        };
    for ifield=1:size(table_map,1);
        if isfield(Filters,table_map{ifield,1})
            aa=table_map{ifield,1};
            bb=table_map{ifield,2};
            cc=table_map{ifield,3};
            for st=1:h
                % in terms of shocks we only save the first-step forecast
                [Filters.(aa){st},Filters.(bb){st},Filters.(cc){st}]=...
                    utils.filtering.squasher(Filters.(aa){st},Filters.(bb){st},init.m_orig);
            end
        end
    end
    
    Fields={'a','att','atT','eta','eta_tt','eta_tlag','epsilon','PAI','PAItt','PAItT';
        'filtered_variables','updated_variables','smoothed_variables',...
        'smoothed_shocks','updated_shocks','filtered_shocks',...
        'smoothed_measurement_errors','filtered_regime_probabilities','updated_regime_probabilities',...
        'smoothed_regime_probabilities'};
    obj.filtering=struct();
    iov=obj.inv_order_var;
    for ifield=1:size(Fields,2)
        main_field=Fields{1,ifield};
        alias=Fields{2,ifield};
        if isfield(Filters,main_field)
            if any(strcmp(main_field,table_map(:,1)))
                % re-order endogenous variables alphabetically
                %-----------------------------------------------
                for reg=1:h
                    Filters.(main_field){reg}=Filters.(main_field){reg}(iov,:,:,:);
                end
            end
            obj.filtering.(alias)=Filters.(main_field);
            if any(strcmp(alias,{'smoothed_variables','smoothed_shocks','smoothed_measurement_errors'}))
                obj.filtering.(['Expected_',alias])=expectation(Filters.PAItT,obj.filtering.(alias));
            elseif any(strcmp(alias,{'updated_variables','updated_shocks'}))
                obj.filtering.(['Expected_',alias])=expectation(Filters.PAItt,obj.filtering.(alias));
            elseif any(strcmp(alias,{'filtered_variables','filtered_shocks'}))
                obj.filtering.(['Expected_',alias])=expectation(Filters.PAI,obj.filtering.(alias));
            end
        end
    end
    %=====================================
    obj=store_probabilities(obj);
    if ~obj.estimation_under_way
        obj=save_filters(obj);
    end
    %=====================================
end

if isempty(LogLik)||isnan(LogLik)||retcode
    % for minimization
    LogLik=-obj.options.estim_penalty;
end

end

function z=recover_positions(obs_id,data_structure,no_more_missing_t,...
    last_good_conditional_observation)

z=@do_it;

    function [ny,occur,obsOccur,no_further_miss,lgco]=do_it(t)
        occur=data_structure(:,t);
        ny=sum(occur); % number of observables to be used in likelihood computation
        obsOccur=obs_id(occur);
        no_further_miss=t>=no_more_missing_t;
        lgco=last_good_conditional_observation(t);
    end
end

function ff=do_one_step_forecast(T,ss,compl,cond_shocks_id,xloc,sig,...
    horizon,is_det_shock)
% y1: forecast
% is_active_shock : location of shocks required to satisfy the constraints.
% sig: perturbation parameter
order=size(T,1);
nstoch=sum(~is_det_shock);
exo_nbr=numel(is_det_shock);
nx=numel(xloc);

ff=@my_one_step;

    function varargout=my_one_step(rt,y0,stoch_shocks,det_shocks)
        % [y1,is_active_shock,retcode,shocks]=my_one_step(rt,y0,shocks)
        if nargin<4
            det_shocks=[];
        end
        y0=struct('y',y0);
        shocks=zeros(exo_nbr,horizon);
        shocks(~is_det_shock,:)=reshape(stoch_shocks,nstoch,horizon);
        if ~isempty(det_shocks)
            span_det=size(det_shocks,2);
            shocks(is_det_shock,1:span_det)=det_shocks;
        end
        nout=nargout;
        if isempty(compl)
            if all(shocks==0) && order==1
                % quick exit if possible
                %------------------------
                y_yss=y0.y-ss{rt};
                y1.y=ss{rt}+T{1,rt}(:,1:nx)*y_yss(xloc)+T{1,rt}(:,nx+1)*sig;
            else
                % if shocks and/or higher orders, do it the hard way
                %----------------------------------------------------
                y1=utils.forecast.one_step_engine(T(:,rt),y0,ss{rt},xloc,sig,...
                    shocks,order);
            end
            outputs={y1,false(1,exo_nbr),0,shocks};
            varargout=cell(1,nout);
            for iarg=1:nout
                varargout{iarg}=outputs{iarg};
            end
        else
            % if restrictions, this is unavoidable
            %--------------------------------------
            [varargout{1:nout}]=utils.forecast.one_step_fbs(T(:,rt),y0,...
                ss{rt},xloc,sig,shocks,order,compl,cond_shocks_id);
        end
        varargout{1}=varargout{1}.y;
    end
end

function E=expectation(probs,vals)
E=0;
for istate=1:numel(vals)
    % E=E+squeeze(bsxfun(@times,permute(probs(istate,:),[3,1,2]),vals{istate}));
    E=E+bsxfun(@times,permute(probs(istate,:),[3,1,2]),vals{istate});
end
% E=squeeze(sum(bsxfun(@times,permute(probs,[3,1,2]),vals),2));
end

function init=re_inflator(init,kf_filtering_level)
m_orig=size(init.a{1},1);
if kf_filtering_level
    is_det_shock=init.is_det_shock;
    exo_nbr=numel(is_det_shock);
    horizon=init.horizon;
    nshocks=exo_nbr*horizon;
    h=numel(init.PAI00);
    nstates=numel(init.state_vars_location);
    for st=1:h
        init.a{st}=[init.a{st};zeros(nshocks,1)];
        init.steady_state{st}=[init.steady_state{st};zeros(nshocks,1)];
        init.P{st}=[init.P{st},zeros(m_orig,nshocks)
            zeros(nshocks,m_orig),eye(nshocks)];
        for io=1:size(init.T,1)
            ncols=size(init.T{io,st},2);
            batch=zeros(nshocks,ncols);
            if io==1
                % place the identities
                batch(:,nstates+2:end)=eye(nshocks);
            end
            init.T{io,st}=[init.T{io,st};batch];
            if io==1
                % collect the Te and Te_det instead of recomputing them
                %-------------------------------------------------------
                Te_all=reshape(init.T{io,st}(:,nstates+2:end),...
                    [m_orig+nshocks,exo_nbr,horizon]);
                init.Te{st}=Te_all(:,~is_det_shock,:);
                init.Te_det{st}=Te_all(:,is_det_shock,:);
                
            end
        end
        % decoupled first order
        %----------------------
        init.Tx{st}=[init.Tx{st};zeros(nshocks,nstates)];
        init.Tsig{st}=[init.Tsig{st};zeros(nshocks,1)];
        
        Te_Te=init.Te{st}(:,:)*init.Te{st}(:,:).';
        Te_Te(1:m_orig,1:m_orig)=init.Te_Te_prime{st};
        init.Te_Te_prime{st}=Te_Te;
    end
end
init.m_orig=m_orig;
init.m=size(init.a{1},1);

end