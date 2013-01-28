function [obj,retcode]=evaluate(obj,varargin)

if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    obj=struct('evaluate_islooping',false,'evaluate_params',[]);
    return
end

% % re-initialize the object
% template=obj.template_obj;
% obj=template;
% obj.template_obj=template;
% take a copy. Should evaluation fail, return the old copy
clean_obj=obj;
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

obj=assign_estimates(obj,params);
% collect the parameters for all regimes
% should probably allow for multiple objects here as well? no
if obj.estimation_under_way
    loc=strcmp('startval',obj.parameters(1,:));
    M=obj.parameters{2,loc};
    loc=strcmp('is_switching',obj.parameters(1,:));
    switching=find(obj.parameters{2,loc});
else
    M=vertcat(obj.parameters.startval);
    switching=find([obj.parameters.is_switching]);
end

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
ss_and_bgp_start_vals=zeros(2*endo_nbr,obj.NumberOfRegimes);

vectorized_code=endo_nbr<=obj.options.vectorized_code_threshold;
dynamic_params=obj.func_handles.dynamic_params;
vectorized_dynamic=obj.func_handles.vectorized_dynamic;
% dynamic=obj.func_handles.dynamic; % both numerical and automatic
% derivatives use the vectorized form of this function, irrespective of the
% flag vectorized_code, which I have to discard, but I am hesitating
endo_exo_derivatives=obj.func_handles.endo_exo_derivatives;
param_derivatives=obj.func_handles.param_derivatives;
planner=obj.func_handles.planner;
vectorized_dynamic_params=obj.func_handles.vectorized_dynamic_params;
definitions=obj.func_handles.definitions;
% steady state functions
ssfunc=obj.func_handles.steady_state_model;
func_ss=obj.func_handles.static;
balanced_growth=obj.func_handles.balanced_growth;
func_jac=obj.func_handles.static_model_derivatives;
static_bgp_model_derivatives=obj.func_handles.static_bgp_model_derivatives;
if obj.options.use_steady_state_model && ~isempty(ssfunc)
    optimopt=optimset('display','none','maxiter',400);
    [ss,retcode]=steady_state_evaluation(obj.parameters,obj.is_unique_steady_state,optimopt,ssfunc);
    if retcode
        retcode=1; % flag on the steady state
        obj=clean_obj;
        if obj(1).options.debug
            decipher_error(retcode)
        end
        return
    end
    % now put the content of ss into ss_and_bgp_start_vals
    ss_and_bgp_start_vals(1:endo_nbr,:)=ss;
    % exogenous auxiliary variables have steady state zero
end
ss_and_bgp_final_vals=ss_and_bgp_start_vals;

derivative_type=obj.options.derivatives;
% steady state for the exogenous
x_ss=zeros(obj.NumberOfExogenous,1);
Lead_lag_incidence=obj.Lead_lag_incidence;
% indices for lags,currents and leads expressed in current period...
indices=[find(Lead_lag_incidence(:,1)>0);(1:endo_nbr)';find(Lead_lag_incidence(:,3)>0)];
%                      leads              current          lags
% threshold for using vectorized code in numerical approximation of the
% jacobian

resid_func=@static_or_bgp_residuals;
last_item=endo_nbr;
if ~isempty(obj.is_stationary_model)
    if ~obj.is_stationary_model
        last_item=2*endo_nbr;
        func_ss=obj.func_handles.balanced_growth;
        func_jac=static_bgp_model_derivatives;
    end
end
NumberOfRegimes=obj.NumberOfRegimes;
def=nan(numel(obj.definitions),NumberOfRegimes);
for ii=NumberOfRegimes:-1:1 % backward to optimize speed
    pp_i=M(:,ii);
    def_i=online_function_evaluator(definitions,pp_i); %#ok<*EVLC>
    if ~obj.is_optimal_policy_model && ~obj.is_imposed_steady_state
        % if the initial guess solves the steady state then proceed. Else
        % try and improve the initial guess through fsolve.
        if obj.is_unique_steady_state
            pp_i=mean(M,2);
        end
        ss_and_bgp_i=ss_and_bgp_start_vals(:,ii);
        if isempty(obj.is_stationary_model)
            % determine stationarity
            
            [is_stationary,ss_and_bgp_i,retcode]=determine_stationarity_status();
            if ~retcode
                obj.is_stationary_model=is_stationary;
                % based on that, set the functions...
                if ~obj.is_stationary_model
                    last_item=2*endo_nbr;
                end
            end
        else
            % solve the steady state
            [ss_and_bgp_i(1:last_item),retcode]=...
                solve_steady_state(ss_and_bgp_i(1:last_item),x_ss,pp_i,def_i,...
                resid_func,obj.is_linear_model,obj.options.optimset);
        end
        if ~retcode
            ss_and_bgp_final_vals(:,ii)=ss_and_bgp_i;
        else
            obj=clean_obj;
            if obj(1).options.debug
                decipher_error(retcode)
            end
            return
        end
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
            if obj(1).options.debug
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
    
    def(:,ii)=def_i;
end
% This steady state is hidden from the users since it might be adjusted or
% expanded when actually solving the model (loose commitment, sticky
% information). It is only when the model is solved that the steady state
% will be pushed into the varendo...
obj.steady_state_and_balanced_growth_path=ss_and_bgp_final_vals;

if ~isempty(def)
    [nrows,ncols]=size(def);
    tmp=mat2cell(def,ones(nrows,1),ncols);
    [obj.definitions(1:nrows).value]=(tmp{:});
end

if ~obj.estimation_under_way && ~obj.is_optimal_policy_model && ...
        ~obj.is_sticky_information_model
    % I wanted to push the steady state right way but sometimes the steady
    % state is not known before the model is totally solved as in the case
    % of optimal policy... and in particular, when the number of endogenous
    % variables changes after the solution
    tmp=mat2cell(obj.steady_state_and_balanced_growth_path,ones(2*endo_nbr,1),NumberOfRegimes);
    tmp0=tmp(1:endo_nbr,:);
    % this thing below works well... even after checking many
    % times :-)
    [obj.varendo(:).det_steady_state]=(tmp0{:});
    tmp0=tmp(endo_nbr+1:end,:);
    [obj.varendo(:).balanced_growth]=(tmp0{:});
end

    function [is_stationary,ys,retcode]=determine_stationarity_status()
        
        is_stationary=[];
        % try stationarity
        
        ys=ss_and_bgp_i;
        [ys(1:endo_nbr),retcode]=solve_steady_state(ss_and_bgp_i(1:endo_nbr),x_ss,pp_i,def_i,...
            @static_or_bgp_residuals,obj.is_linear_model,obj.options.optimset);
        
        if retcode
            % try nonstationarity
            func_ss=balanced_growth;
            func_jac=static_bgp_model_derivatives;
            [ys,retcode]=solve_steady_state(ss_and_bgp_i,x_ss,pp_i,def_i,...
                @static_or_bgp_residuals,obj.is_linear_model,obj.options.optimset);
            if ~retcode
                is_stationary=false;
                if ~obj.options.evaluate_islooping
                    disp([mfilename,':: model is not mean-stationary but allows for a BALANCED GROWTH PATH'])
                end
            else
                if ~obj.options.evaluate_islooping
                    disp([mfilename,':: stationarity status could not be determined'])
                end
            end
        else
            is_stationary=true;
            if ~obj.options.evaluate_islooping
                disp([mfilename,':: model found to be MEAN-stationary'])
            end
        end
    end

    function [lhs,Jac,retcode]=static_or_bgp_residuals(y_,x_,param_,def_)
        % if the structural model is defined in terms of steady states, then the
        % the word ss will appear in the static model. But then, those values are
        % the same as those found in y.
        % In some cases, the model is nonstationary so that only the BGP is
        % solvable in this case, the y vector also contains the balanced-growth
        % values. copying the whole vector to ss does not harm since only the
        % relevant elements in the upper part of y will be picked up.
        
        retcode=0;
        
        ss_=y_;
        lhs=online_function_evaluator(func_ss,y_,x_,param_,ss_,def_);
%         if any(any(isnan(lhs)))||(any(any(isnan(lhs))))
%             retcode=?
%             if retcode, Jac=[];
%             end
%         end
        if nargout>1
            Jac=online_function_evaluator(func_jac,y_,x_,param_,ss_,def_);
%             if any(any(isnan(Jac)))||(any(any(isnan(Jac))))
%                 retcode=?
%             end
        end
    end
end

function [ss,retcode]=steady_state_evaluation(param_obj,unique_ss,options,ssfunc) %#ok<INUSL>

if nargin<3
    options=optimset('display','none','maxiter',400,'tolfun',1e-6); %#ok<NASGU>
    if nargin<2
        unique_ss=true;
    end
end

retcode=0;
if isa(param_obj,'rise_param')
    par_vals=vertcat(param_obj.startval);
else
    loc= strcmp('startval',param_obj(1,:));
    par_vals=vertcat(param_obj{2,loc});
end
number_of_regimes=size(par_vals,2);
if unique_ss
    par_vals=mean(par_vals,2);
end
ss=[];
for ii=1:number_of_regimes
    if ii==1 || ~unique_ss
        param=par_vals(:,ii); 
        % in ssfunc, the definitions have been substituted already in
        % load_functions and so all is needed here is the parameters
        y=online_function_evaluator(ssfunc,[],[],param,[],[]);
    end
    if retcode
        return %#ok<UNRCH>
    end
    if ii==1
        ss=nan(numel(y),number_of_regimes);
    end
    ss(:,ii)=y; %#ok<AGROW>
end

end

