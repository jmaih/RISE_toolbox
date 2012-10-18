function [obj,LogLik,Incr,retcode]=filter(obj,varargin)

if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    obj=markov_switching_kalman_filter();
    return
end

params=[];
if ~isempty(varargin)
	tmp=find(strcmp('evaluate_params',varargin));
	if ~isempty(tmp)
		params=varargin{tmp+1};
		varargin(tmp+[0,1])=[];
	end
	obj=set_options(obj,varargin{:});
end

if ~obj.data_are_loaded
    obj=obj.load_data;
end

Incr=[];

resolve_flag=isempty(obj.T)||~isempty(params);

% solve the object
if resolve_flag
    [obj,retcode]=obj.solve('evaluate_params',params);
    if retcode
        LogLik=-obj.options.Penalty;
		if obj(1).options.debug
			decipher_error(retcode)
		end
        return
    end
end

n_endo=obj.NumberOfEndogenous(2);
% extract the steady state from the ss_and_bgp
SS=obj.steady_state_and_balanced_growth_path(1:n_endo,:);

% extract the state and transition matrices
T=obj.T;
R=obj.R;
Q={obj.Q,[],[]};
if  obj.is_endogenous_switching_model
    % take a handle to a private function. we can't access it otherwise
    Q{2}=obj.func_handles.transition_matrix;
%    Q{2}=@transition_matrix_evaluation;
%    transition_script=obj.func_handles.transition_matrix;

    % collect the parameters for all regimes
%     M=vertcat(obj.parameters.startval);
    M=obj.parameters{2,end};
    % transition matrix should be invariant. And so, hopefully, the first
    % argument to Q{2} could be the updated data in any state.
    % order of the input arguments is y,x,ss,param
    Q{3}={[],SS(:,1),mean(M,2)}; % remaining arguments of Q{2} after the first one
%    Q{3}={ss(:,1),mean(M,2),transition_script}; % remaining arguments of Q{2}
    % [obj.Q,retcode]=transition_matrix_evaluation(ss_i,ss_i,M(:,1),transition_matrix); % should be invariant
end

H=obj.H;% measurement errors

% observations
data=vertcat(obj.varobs.value);
% location of observables in the state
data_id=vertcat(obj.varobs.id);
% conditional observations
cond_data=vertcat(obj.cond_varobs.value);
% Now the subtelty: if order is 0, then there is no future information
% even if we later on I go back to estimating the expansion order, I
% won't have to change anything here. I just need to assign to order
% its estimated value... but then a parameter with a specific name
% will have to be defined and we will have to impose that the parameter
% is constant, even if in theory, one could well assume that in a stable
% regime, people have a better foresight of the future than in a volatile
% regime...
if obj.options.order==0
    cond_data=nan*cond_data;
end

% location of observables in the state
cond_data_id=vertcat(obj.cond_varobs.id); % probably cat(1,... would do too
% if there is no conditional information, the expansion order should be
% exactly 1
if isempty(cond_data_id) && size(R,3)~=1
    error([mfilename,':: with no conditional information, the expansion order (currently ',int2str(size(R,3)),') should be 1'])
end
% evaluate the likelihood

OMG=[];
Hypothesis=obj.options.forecast_conditional_hypothesis;

%% Check the waters to detect whether we need to transform the system prior to filtering
% There will be as many preconditioning matrices as the number of regimes.
% one could also simply do as Tao and Dan, that is to use a diffuse
% initialization with high variance and presample no matter whether the
% variables are stationary or not 

%{
precondition=false(1,obj.NumberOfRegimes);
Us=cell(1,obj.NumberOfRegimes);
P0=repmat(diag(lyapunov_diffuse_factor*ones(n_endo,1)),[1,1,obj.NumberOfRegimes]);
for ii=1:obj.NumberOfRegimes
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

%%
% Partition R. 
% list of deterministic shocks
det_shocks={obj.varobs_exo.name};
% locate those shocks among the exogenous
det_shocks_id=locate_variables(det_shocks,{obj.varexo.name});
is_deterministic=false(1,obj.NumberOfExogenous);
is_deterministic(det_shocks_id)=true;
R_det=reshape(R(:,det_shocks_id,1,:),[n_endo,numel(det_shocks_id),obj.NumberOfRegimes]); %squeeze() take only the first order
% constant term reshaped to match R_det
risk=reshape(vertcat(obj.varendo.risk),[n_endo,1,obj.NumberOfRegimes]);
R_det=cat(2,risk,R_det);
if all(all(all(R_det==0)))
    R_det=[];
end
% Finally the R for the truly exogenous variables
R=R(:,~is_deterministic,:,:);
if ~isempty(R_det)
    % load deterministic shock information
    exo_data=vertcat(obj.varobs_exo.value);
    % add a 1 as the first row of the deterministic data
    exo_data=[ones(1,size(data,2));exo_data];
    WB={R_det,exo_data};
else
    WB={};
end

if obj.options.kf_filtering_level || obj.is_endogenous_switching_model
    if isempty(cond_data_id)
        [LogLik,Incr,retcode,Filters] = markov_switching_kalman_filter(data_id,data,...
            T,R,SS,Q,H,obj.options,WB);
    else
        [LogLik,Incr,retcode,Filters] = markov_switching_kalman_filter_real_time(data_id,data,...
            cond_data_id,cond_data,OMG,Hypothesis,...
            T,R,SS,Q,H,obj.options,WB);
    end
    if ~retcode
        Fields={'a','att','alphat','eta','epsilon','BIGPAI','BIGPAI_tt','BIGPAI_tT';
            'filtered_variables','updated_variables','smoothed_variables','smoothed_shocks',...
            'smoothed_measurement_errors','filtered_probabilities','updated_probabilities',...
            'smoothed_probabilities'};
        obj.Filters=struct();
        for ifield=1:size(Fields,2)
            if isfield(Filters,Fields{1,ifield})
                obj.Filters.(Fields{2,ifield})=Filters.(Fields{1,ifield});
                if ismember(Fields{2,ifield},{'smoothed_variables','smoothed_shocks','smoothed_measurement_errors'})
                    obj.Filters.(['Expected_',Fields{2,ifield}])=expectation(Filters.BIGPAI_tT,obj.Filters.(Fields{2,ifield}));
                elseif strcmp(Fields{2,ifield},'updated_variables')
                    obj.Filters.(['Expected_',Fields{2,ifield}])=expectation(Filters.BIGPAI_tt,obj.Filters.(Fields{2,ifield}));
                elseif strcmp(Fields{2,ifield},'filtered_variables')
                    obj.Filters.(['Expected_',Fields{2,ifield}])=expectation(Filters.BIGPAI,obj.Filters.(Fields{2,ifield}));
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
    % N.B: the endogenous probability script is written as a function of
    % the grand state (i.e. all the endogenous variables) and so, it is
    % dangerous to reduce the state in that case. The drawback is that
    % likelihood computation becomes slower. Amending this would amount to
    % having a separate script for the transition matrix for the specific
    % case where we don't need all the smoothed variables. This will have
    % to be done from the parser, when reading the list of the observables.
    % Or something along those lines...
else
    [state_cols,new_data_id,new_cond_data_id]=minimum_state_for_estimation(T,data_id,cond_data_id);
    if ~isempty(R_det)
        WB{1}=WB{1}(state_cols,:,:);
    end
    if isempty(cond_data_id)
        [LogLik,Incr,retcode] = markov_switching_kalman_filter(new_data_id,data,T(state_cols,state_cols,:),...
            R(state_cols,:,:,:),SS(state_cols,:),Q,H,obj.options,WB);
    else
        [LogLik,Incr,retcode] = markov_switching_kalman_filter_real_time(new_data_id,data,...
            new_cond_data_id,cond_data,OMG,Hypothesis,T(state_cols,state_cols,:),...
            R(state_cols,:,:,:),SS(state_cols,:),Q,H,obj.options,WB);
    end
end
if isnan(LogLik)||retcode
    % for maximization
    LogLik=-obj.options.Penalty;
end

function E=expectation(probs,vals)
E=squeeze(sum(bsxfun(@times,permute(probs,[3,1,2]),vals),2));

function varargout=minimum_state_for_estimation(state,varargin)
% this function reduces the state size to accelerate estimation
% it returns the indices of the union of state variables and observable
% ones
if ~islogical(state)
    tmp=state;
    state=any(tmp(:,:,1));
    for ii=2:size(tmp,3)
        state=state | any(tmp(:,:,ii));
    end
end
n=numel(state);
newstate=state;
for ii=1:length(varargin)
    newstate(varargin{ii})=true;
end
varargout{1}=newstate;
for ii=1:length(varargin)
    newobs=false(1,n);
    newobs(varargin{ii})=true;
    newobs(~newstate)=[];
    newobs=find(newobs);
    [~,tags]=sort(varargin{ii});
    newobs(tags)=newobs;
    varargout{ii+1}=newobs;
end
