function [obj,retcode]=evaluate(obj,varargin)

if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    obj=struct('evaluate_islooping',false,'evaluate_params',[]);
    return
end

% assign the new estimates
params=[];
if ~isempty(varargin)
    tmp=find(strcmp('evaluate_params',varargin));
    if ~isempty(tmp)
        params=varargin{tmp+1};
        varargin(tmp+[0,1])=[];
    end
    obj=set_options(obj,varargin{:});
end

% take a copy. Should evaluation fail, return the old copy
clean_obj=obj;

obj=assign_estimates(obj,params);
% collect the parameters for all regimes
% should probably allow for multiple objects here as well? no

% compute the steady state and possibly update the parameters if there is a
% steady state file in which parameters are functions of variables instead
% of the other way around
[obj,ss_and_bgp_final_vals,retcode]=compute_steady_state(obj);
if retcode
    obj=clean_obj;
    if obj(1).options.debug
        decipher_error(retcode)
    end
    return
end

loc=strcmp('startval',obj.parameters_image(1,:));
M=obj.parameters_image{2,loc};
loc=strcmp('is_switching',obj.parameters_image(1,:));
switching=find(obj.parameters_image{2,loc});

if ~obj.estimation_under_way
    parameter_restrictions=obj.func_handles.parameter_restrictions;
    retcode=parameter_restrictions(M);
    if retcode
        if obj(1).options.debug
            decipher_error(retcode)
        end
        return
    end
end

% measurement errors
measure_flag=false;
if ~isempty(obj.measurement_errors_restrictions)
    obj.H=zeros(obj.NumberOfObservables(1),obj.NumberOfObservables(1),obj.NumberOfRegimes);
    measure_flag=true;
    Restrictions=obj.measurement_errors_restrictions;
end

% we compute the steady state before solving the model...
% default a steady state file, the steady_state solver below, the initial
% guess is zero
endo_nbr=obj.NumberOfEndogenous(1);

% % vectorized_code=endo_nbr<=obj.options.vectorized_code_threshold;
% % dynamic_params=obj.func_handles.dynamic_params;
vectorized_dynamic=obj.func_handles.vectorized_dynamic;
% dynamic=obj.func_handles.dynamic; % both numerical and automatic
% derivatives use the vectorized form of this function, irrespective of the
% flag vectorized_code, which I have to discard, but I am hesitating
endo_exo_derivatives=obj.func_handles.endo_exo_derivatives;
param_derivatives=obj.func_handles.param_derivatives;
planner=obj.func_handles.planner;
vectorized_dynamic_params=obj.func_handles.vectorized_dynamic_params;

derivative_type=obj.options.derivatives;
% steady state for the exogenous
x_ss=zeros(obj.NumberOfExogenous,1);
Lead_lag_incidence=obj.Lead_lag_incidence;
% indices for lags,currents and leads expressed in current period...
indices=[find(Lead_lag_incidence(:,1)>0);(1:endo_nbr)';find(Lead_lag_incidence(:,3)>0)];
%                      leads              current          lags
% threshold for using vectorized code in numerical approximation of the
% jacobian

NumberOfRegimes=obj.NumberOfRegimes;
% get the value of the definitions. These were evaluated in compute_steady_state above
def=vertcat(obj.definitions.value);
def_i=[];
for ii=NumberOfRegimes:-1:1 % backward to optimize speed
    if ~isempty(def)
        def_i=def(:,ii); %#ok<*EVLC>
    end
    
    % with the steady state in hand, we can go ahead and compute the
    % jacobian conditional on the steady state.
    % first, get the steady state for all the variables in the incidence
    % matrix
    ss_i=ss_and_bgp_final_vals(1:endo_nbr,ii);
    % build the jacobian
    
    switch derivative_type
        case 'symbolic'
            J=online_function_evaluator(endo_exo_derivatives,ss_i(indices),x_ss,M(:,ii),ss_i,def_i);
            if ~isempty(switching)
                JP=online_function_evaluator(param_derivatives,ss_i(indices),x_ss,M(:,ii),ss_i,def_i); %#ok<*GTARG>
                J=[J,JP]; %#ok<AGROW>
            end
        case 'numerical'
            z=[ss_i(indices);x_ss];
            myfun=@(xx)online_function_evaluator(vectorized_dynamic,xx,M(:,ii),ss_i,def_i);
            J=rise_numjac(myfun,z);
            if ~isempty(switching)
                % now the z vector has to be a matrix in order to avoid a
                % breakdown occurring when no parameter is present in some
                % equations.
                z=z(:,ones(1,size(M,1)));
                myfun=@(xx)online_function_evaluator(vectorized_dynamic_params,xx,z,ss_i);
                JP=rise_numjac(myfun,M(:,ii));
                J=[J,JP]; %#ok<AGROW>
            end
        case 'automatic'
            z=[ss_i(indices);x_ss];
            zz=rise_nad(z);
            J=online_function_evaluator(vectorized_dynamic,zz,M(:,ii),ss_i,def_i);
            if ~isempty(switching)
                JP=online_function_evaluator(vectorized_dynamic_params,rise_nad(M(:,ii)),z,ss_i);
                J=[J,JP]; %#ok<AGROW>
            end
            J=full(J);
    end
    if any(any(isnan(J)))
        retcode=2; % nans in jacobian
        obj=clean_obj;
        if obj(1).options.debug
            decipher_error(retcode)
        end
        return
    end
    
    if ii==NumberOfRegimes
        % A nice result is that at first-order approximation, the jacobians
        % are just reweighted by the probabilities evaluated at the steady
        % state.
        [obj.Q,retcode]=online_function_evaluator(obj.func_handles.transition_matrix,ss_i,[],M(:,1),ss_i,def_i);
        if retcode
            obj=clean_obj;
            if obj.options.debug
                decipher_error(retcode)
            end
            return
        end
        % initialize
        Aminus_=zeros(obj.NumberOfEquations,endo_nbr);
        A0_=Aminus_;
        Aplus_=Aminus_;
        % lags
        [aa_,~,cc_]=find(obj.Lead_lag_incidence(:,3));
        % current
        [aa0,~,cc0]=find(obj.Lead_lag_incidence(:,2));
        % leads
        [aap,~,ccp]=find(obj.Lead_lag_incidence(:,1));
        % shocks
        cc_shocks=nnz(obj.Lead_lag_incidence)+(1:obj.NumberOfExogenous);
        % switching parameters
        cc_switch=nnz(obj.Lead_lag_incidence)+obj.NumberOfExogenous+switching;
        obj.C=nan(obj.NumberOfEquations,numel(switching),NumberOfRegimes);
    end
    Aminus_(:,aa_)=J(:,cc_);
    A0_(:,aa0)=J(:,cc0);
    Aplus_(:,aap)=J(:,ccp);
    [obj.Aminus(:,:,ii),obj.A0(:,:,ii),obj.Aplus(:,:,ii)]=deal(Aminus_,A0_,Aplus_);
    
    obj.B(:,:,ii)=J(:,cc_shocks);
    
    % switching parameters
    obj.C(:,:,ii)=J(:,cc_switch);
    
    % loss function weights
    if obj.is_optimal_policy_model || obj.is_optimal_simple_rule_model
        [obj.planner_loss(1,ii),...
            obj.planner_commitment(1,ii),...
            obj.planner_discount(1,ii),W1]=online_function_evaluator(...
            planner,ss_i,x_ss,M(:,ii),ss_i,def_i);
        if obj.options.debug
            plan_fd=@(zz)online_function_evaluator(planner,zz,x_ss,M(:,ii),ss_i,def_i);
            W2=finite_difference_hessian(plan_fd,ss_i);
            disp('discrepancies in planner computations')
            disp(max(max(abs(W1-W2))))
        else
            if strcmp(derivative_type,'symbolic')
                obj.W(:,:,ii)=W1;
            else
                plan_fd=@(zz)online_function_evaluator(planner,zz,x_ss,M(:,ii),ss_i,def_i);
                W2=finite_difference_hessian(plan_fd,ss_i);
                obj.W(:,:,ii)=W2;
            end
        end
    end
    if measure_flag
        tmp=zeros(obj.NumberOfObservables(1),1);
        tmp(Restrictions(:,1))=M(Restrictions(:,2),ii).^2;
        obj.H(:,:,ii)=diag(tmp);
    end
end

end


