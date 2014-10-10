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
% Description
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
    obj=msre_linear_filter();%markov_switching_kalman_filter();
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

if ~isempty(varargin)
	obj=set(obj,varargin{:});
end

if ~obj.data_are_loaded
    obj=obj.load_data;
end

resolve_flag=isempty(obj.solution)||isempty(obj.solution.Tz{1});
%
%% solve the object
if resolve_flag
    [obj,retcode]=obj.solve();
    if retcode
        LogLik=-obj.options.estim_penalty;
		if obj(1).options.debug
			utils.error.decipher(retcode)
		end
        return
    end
end
h=obj.markov_chains.regimes_number;

% extract the steady state
SS=obj.solution.ss;

ov=obj.order_var.after_solve;
iov=obj.inv_order_var.after_solve;
zpb=obj.locations.after_solve.z.pb;
tpb=obj.locations.after_solve.t.pb;
ze_0=obj.locations.after_solve.z.e_0;
zsig=obj.locations.after_solve.z.sig;
endo_nbr=obj.endogenous.number(end);

% extract the state and transition matrices
T=cell(1,h); 
R=T;
% risk term
%----------
risk=T;
reorder_rows=false;
Tall=zeros(endo_nbr);
for ireg=1:h
    tmp=obj.solution.Tz{ireg};
    if reorder_rows
        % re-order the steady state
        SS{ireg}=SS{ireg}(ov);
        % re-order the transition and shock matrices
        tmp=tmp(ov,:);
    end
    Tall(:,tpb)=tmp(:,zpb);
    T{ireg}=Tall;
    if ~reorder_rows
        % reorder the columns then
        T{ireg}=T{ireg}(:,iov);
    end
    R{ireg}=tmp(:,ze_0);
    risk{ireg}=obj.options.simul_sig*tmp(:,zsig);
end
clear Tall
% extract data and put them in the order_var order as well
%---------------------------------------------------------
data=obj.data;
if reorder_rows
    data.varobs_id=iov(data.varobs_id);
end

Qfunc=prepare_transition_routine(obj);

H=obj.solution.H;% measurement errors

% deterministic shocks
det_shocks=obj.exogenous.is_observed;

% load deterministic shock data
%------------------------------
exo_data=[ones(1,size(obj.data.x,2));obj.data.x];

% constant and deterministic terms
%---------------------------------
I=speye(endo_nbr);
State_trend=cell(1,h);
for istate=1:h
    const_st=SS{istate};
    if any(const_st)
        const_st=(I-T{istate})*const_st;
    end
    const_st=const_st+risk{istate};
    Coef_det_st=[const_st,R{istate}(:,det_shocks,1)];
    State_trend{istate}=Coef_det_st*exo_data;
    % Trim R for the stochastic exogenous variables
    %----------------------------------------------
    R{istate}(:,det_shocks,:)=[];
end

% initialization
%---------------

syst=struct('T',{T},'R',{R},'H',{H},'Qfunc',{Qfunc});

data_trend=[];
[LogLik,Incr,retcode,Filters]=msre_linear_filter(syst,data,data_trend,State_trend,SS,risk,obj.options);
if  obj.options.kf_filtering_level && ~retcode
    Fields={'a','att','atT','eta','epsilon','PAI','PAItt','PAItT';
        'filtered_variables','updated_variables','smoothed_variables','smoothed_shocks',...
        'smoothed_measurement_errors','filtered_regime_probabilities','updated_regime_probabilities',...
        'smoothed_regime_probabilities'};
    obj.filtering=struct();
    if reorder_rows
        % things need to be re-ordered here
        keyboard
    end
    for ifield=1:size(Fields,2)
        if isfield(Filters,Fields{1,ifield})
            obj.filtering.(Fields{2,ifield})=Filters.(Fields{1,ifield});
            if ismember(Fields{2,ifield},{'smoothed_variables','smoothed_shocks','smoothed_measurement_errors'})
                obj.filtering.(['Expected_',Fields{2,ifield}])=expectation(Filters.PAItT,obj.filtering.(Fields{2,ifield}));
            elseif strcmp(Fields{2,ifield},'updated_variables')
                obj.filtering.(['Expected_',Fields{2,ifield}])=expectation(Filters.PAItt,obj.filtering.(Fields{2,ifield}));
            elseif strcmp(Fields{2,ifield},'filtered_variables')
                obj.filtering.(['Expected_',Fields{2,ifield}])=expectation(Filters.PAI,obj.filtering.(Fields{2,ifield}));
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
if isnan(LogLik)||retcode
    % for minimization
    LogLik=-obj.options.estim_penalty;
end

function E=expectation(probs,vals)
E=0;
for istate=1:numel(vals)
% E=E+squeeze(bsxfun(@times,permute(probs(istate,:),[3,1,2]),vals{istate}));
E=E+bsxfun(@times,permute(probs(istate,:),[3,1,2]),vals{istate});
end
% E=squeeze(sum(bsxfun(@times,permute(probs,[3,1,2]),vals),2));

% function varargout=minimum_state_for_estimation(state,varargin)
% % this function reduces the state size to accelerate estimation
% % it returns the indices of the union of state variables and observable
% % ones
% if ~islogical(state)
%     tmp=state;
%     % The step below is critical for speed, though it may add some noise
%     % when computing the filters for all the endogenous variables
%     tmp(abs(tmp)<1e-9)=0; 
%     state=any(tmp(:,:,1));
%     for ii=2:size(tmp,3)
%         state=state | any(tmp(:,:,ii));
%     end
% end
% n=numel(state);
% n_varg=length(varargin);
% newstate=state;
% for ii=1:n_varg
%     newstate(varargin{ii})=true;
% end
% varargout=cell(1,n_varg+1);
% varargout{1}=newstate;
% for ii=1:n_varg
%     newobs=false(1,n);
%     newobs(varargin{ii})=true;
%     newobs(~newstate)=[];
%     newobs=find(newobs);
%     [~,tags]=sort(varargin{ii});
%     newobs(tags)=newobs;
%     varargout{ii+1}=newobs;
% end
