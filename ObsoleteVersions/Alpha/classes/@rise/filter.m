function [obj,LogLik,Incr,retcode]=filter(obj,varargin)

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
	obj=set_options(obj,varargin{:});
end

if ~obj.data_are_loaded
    obj=obj.load_data;
end

resolve_flag=isempty(obj.solution.m_x{1});
%
%% solve the object
if resolve_flag
    [obj,retcode]=obj.solve();
    if retcode
        LogLik=-obj.options.Penalty;
		if obj(1).options.debug
			decipher_error(retcode)
		end
        return
    end
end

% extract the steady state
SS=obj.solution.ss;

% extract the state and transition matrices
T=obj.solution.m_x;
R=obj.solution.m_e;
Q={obj.solution.Q,[],[]};
% risk term
%----------
risk=obj.solution.m_sig;

if  obj.is_endogenous_switching_model
    % take a handle to a private function. we can't access it otherwise
    Q{2}=obj.func_handles.transition_matrix;
%    Q{2}=@transition_matrix_evaluation;
%    transition_script=obj.func_handles.transition_matrix;

    % collect the parameters for all regimes
%     M=vertcat(obj.parameters.startval);
    M=obj.parameter_values;
    % transition matrix should be invariant. And so, hopefully, the first
    % argument to Q{2} could be the updated data in any state.
    % order of the input arguments is y,x,param,ss
    Q{3}={[],mean(M,2),SS{1}}; % remaining arguments of Q{2} after the first one
%    Q{3}={ss(:,1),mean(M,2),transition_script}; % remaining arguments of Q{2}
    % [obj.Q,retcode]=transition_matrix_evaluation(ss_i,ss_i,M(:,1),transition_matrix); % should be invariant
end

H=obj.solution.H;% measurement errors

% % observations
% cond_data=vertcat(obj.cond_varobs.value);
% % Now the subtelty: if order is 0, then there is no future information
% % even if we later on I go back to estimating the expansion order, I
% % won't have to change anything here. I just need to assign to order
% % its estimated value... but then a parameter with a specific name
% % will have to be defined and we will have to impose that the parameter
% % is constant, even if in theory, one could well assume that in a stable
% % regime, people have a better foresight of the future than in a volatile
% % regime...
% if obj.options.order==0
%     cond_data=nan*cond_data;
% % end
% 
% % location of observables in the state
% cond_data_id=vertcat(obj.cond_varobs.id); % probably cat(1,... would do too
% % if there is no conditional information, the expansion order should be
% % exactly 1
% if isempty(cond_data_id) && size(R,3)~=1
%     error([mfilename,':: with no conditional information, the expansion order (currently ',int2str(size(R,3)),') should be 1'])
% end
% % evaluate the likelihood
% 
% OMG=[];
% Hypothesis=obj.options.forecast_conditional_hypothesis;

%% Check the waters to detect whether we need to transform the system prior to filtering
% There will be as many preconditioning matrices as the number of regimes.
% one could also simply do as Tao and Dan, that is to use a diffuse
% initialization with high variance and presample no matter whether the
% variables are stationary or not 

%{
precondition=false(1,obj.markov_chains.regimes_number);
Us=cell(1,obj.markov_chains.regimes_number);
P0=repmat(diag(lyapunov_diffuse_factor*ones(n_endo,1)),[1,1,obj.markov_chains.regimes_number]);
for ii=1:obj.markov_chains.regimes_number
    if max(abs(eig(T(:,:,ii))))>1-tol_________________
        precondition(ii)=true;
        [Ts,RRs,Us{ii},stable]=lower_triangularize(T(:,:,ii),RR(:,:,ii),obj.options);
        [Vs,retcode]=solve_triangular_lyapunov(Ts(stable,stable),Qs(stable,stable));
        if retcode
            LogLik=-obj.options.Penalty;
			if obj(1).options.debug
				decipher_error(retcode)
			end
            return
        end
        P0(stable,stable)=Vs;
        T(:,:,ii)=Ts;
        RR(:,:,ii)=RRs;
    else
    end
end
%}

h=obj.markov_chains.regimes_number;
% deterministic shocks
det_shocks=obj.exogenous.is_observed;

% load deterministic shock data
%------------------------------
exo_data=[ones(1,size(obj.data.x,2));obj.data.x];

% constant and deterministic terms
%---------------------------------
I=speye(size(T{1},1));
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

syst=struct('T',{T},'R',{R},'H',{H},'Q0',{Q});

data_trend=[];
[LogLik,Incr,retcode,Filters]=msre_linear_filter(syst,obj.data,data_trend,State_trend,SS,risk,obj.options);
if  obj.options.kf_filtering_level && ~retcode
    Fields={'a','att','atT','eta','epsilon','PAI','PAItt','PAItT';
        'filtered_variables','updated_variables','smoothed_variables','smoothed_shocks',...
        'smoothed_measurement_errors','filtered_regime_probabilities','updated_regime_probabilities',...
        'smoothed_regime_probabilities'};
    obj.filtering=struct();
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
    LogLik=-obj.options.Penalty;
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
