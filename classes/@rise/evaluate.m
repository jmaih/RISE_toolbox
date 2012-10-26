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

ssfunc=obj.func_handles.steady_state_model;
vectorized_dynamic=obj.func_handles.vectorized_dynamic;
dynamic=obj.func_handles.dynamic;
vectorized_code=endo_nbr<=obj.options.vectorized_code_threshold;
static=obj.func_handles.static;
balanced_growth=obj.func_handles.balanced_growth;
planner=obj.func_handles.planner;
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

if ~isempty(obj.is_stationary_model)
    if obj.is_stationary_model
        resid_func=@steady_state_residuals;
        last_item=endo_nbr;
        static_func=obj.func_handles.static;
    else
        resid_func=@balanced_growth_path_residuals;
        last_item=2*endo_nbr;
        static_func=obj.func_handles.balanced_growth;
    end
end
NumberOfRegimes=obj.NumberOfRegimes;
for ii=NumberOfRegimes:-1:1 % backward to optimize speed
    if ~obj.is_optimal_policy_model && ~obj.is_imposed_steady_state
        % if the initial guess solves the steady state then proceed. Else
        % try and improve the initial guess through fsolve.
        pp_i=M(:,ii);
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
                if obj.is_stationary_model
                    resid_func=@steady_state_residuals;
                    last_item=endo_nbr;
                    static_func=static;
                else
                    resid_func=@balanced_growth_path_residuals;
                    last_item=2*endo_nbr;
                    static_func=balanced_growth;
                end
            end
        else
            % solve the steady state
            [ss_and_bgp_i(1:last_item),retcode]=...
                solve_steady_state(ss_and_bgp_i(1:last_item),x_ss,pp_i,...
                resid_func,static_func,obj.is_linear_model,obj.options.optimset);
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
 
    J=eval(obj.func_handles.derivatives,ss_i(indices),x_ss,ss_i,M(:,ii));
    J=J{1};
%    z=[ss_i(indices);x_ss;M(switching,ii)];
%    switch derivative_type
%         case 'numerical'
%             J=numerical_jacobian(vectorized_dynamic,z,vectorized_code,ss_i,M(:,ii));
%         case 'automatic'
%             if ii==1
%                 zmad=automatic(z);
%             else
%                 zmad=set(zmad,'x',z(:));
%             end
%             J=get(dynamic(zmad,ss_i,M(:,ii)),'dx');
% %             J=full(first_order_derivatives(dynamic,z,ss_i,M(:,ii)));
%     end
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
        % Ideally the transition matrix should not be a function of something that switches...
        % but then, does that mean the steady state should be constant?
        [obj.Q,retcode]=kron(obj.func_handles.transition_matrix,ss_i,ss_i,[],M(:,1));
%         [obj.Q,retcode]=transition_matrix_evaluation(ss_i,ss_i,M(:,1),...
%             obj.func_handles.transition_matrix);
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
        [~,obj.planner_commitment(1,ii),...
            obj.planner_discount(1,ii),~,W1]=policy_evaluation(ss_i,ss_i,M(:,ii),planner);
        if obj.options.debug
            W2=finite_difference_hessian(@policy_evaluation,ss_i,ss_i,M(:,ii),planner);
            disp('discrepancies in planner computations')
            disp(W1-W2)
        else
            if strcmp(derivative_type,'symbolic')
                obj.W(:,:,ii)=W1;
            else
                obj.W(:,:,ii)=finite_difference_hessian(@policy_evaluation,ss_i,ss_i,M(:,ii),planner);
            end
        end
    end
    if measure_flag
        tmp=zeros(obj.NumberOfObservables(1),1);
        tmp(Restrictions(:,1))=M(Restrictions(:,2),ii).^2;
        obj.H(:,:,ii)=diag(tmp);
    end
end
% This steady state is hidden from the users since it might be adjusted or
% expanded when actually solving the model (loose commitment, sticky
% information). It is only when the model is solved that the steady state
% will be pushed into the varendo...
obj.steady_state_and_balanced_growth_path=ss_and_bgp_final_vals;

%% Sparse everything, unfortunately it does not work on arrays of more than 2 dimensions
%SparseList={'Aminus','A0','Aplus','Q','B','H'};
%for kk=1:numel(SparseList)
%	obj.(SparseList{kk})=sparse(obj.(SparseList{kk}));
%end

    function [is_stationary,ys,retcode]=determine_stationarity_status()
        
        is_stationary=[];
        % try stationarity
        ys=ss_and_bgp_i;
        [ys(1:endo_nbr),retcode]=solve_steady_state(ss_and_bgp_i(1:endo_nbr),x_ss,pp_i,...
            @steady_state_residuals,static,obj.is_linear_model,obj.options.optimset);
        if retcode
            % try nonstationarity
            [ys,retcode]=solve_steady_state(ss_and_bgp_i,x_ss,pp_i,...
                @balanced_growth_path_residuals,balanced_growth,obj.is_linear_model,obj.options.optimset);
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
end

function [objective,commitment,discount,der1,der2]=policy_evaluation(y,ss,param,planner) %#ok<STOUT,INUSL>
eval(planner);
end

function J=numerical_jacobian(objfun,y0,vectorized_code,varargin)
if nargin<3 || isempty(vectorized_code)
    vectorized_code=true;
end

% one evaluation of the objective just to get the number of equations.
% Hopefully this is not expensive. I would have to do this if I had to
% evaluate one-sided derivatives...
f0=objfun(y0,varargin{:});
NumberOfEquations=size(f0,1);

delta_y=sqrt(eps)*max(abs(y0),1);
temp = max(abs(y0),1)+delta_y; delta_y=temp-max(abs(y0),1);
% above we make sure that h and x+h are representable in floating point
% precision so that the difference between x+h and x is exactly h. This is
% accomplished by placing their values into and out of memory as
% follows: h = sqrt(eps) * x, temp = x + h and h = temp ? x.

ysize=numel(y0);
FPLUS=nan(NumberOfEquations,ysize);
FMINUS=nan(NumberOfEquations,ysize);
% we now build each column of the jacobian using numerical derivatives
if vectorized_code
    ones_y=ones(1,ysize);
    y=y0(:,ones_y)+diag(delta_y);
    FPLUS=objfun(y,varargin{:});
    y=y0(:,ones_y)-diag(delta_y);
    FMINUS=objfun(y,varargin{:});
else
    % direction=1:xsize+ysize;
    direction=ysize:-1:1;
    for ii=direction % evaluate backward for faster allocation
        y=y0(:,ones(2,1));
        y(ii,:)=y(ii,:)+[delta_y(ii),-delta_y(ii)];
        FPLUS(:,ii)=objfun(y(:,1),varargin{:});
        FMINUS(:,ii)=objfun(y(:,2),varargin{:});
    end
end
% central differences more accurate...
delta=2*delta_y';
J=(FPLUS-FMINUS)./delta(ones(NumberOfEquations,1),:);
end

function res=derivatives_evaluation(ss,x,param,indices,deriv_func) %#ok<STOUT,INUSL>
y=ss(indices);
y=y(:); %#ok<NASGU>
eval(deriv_func)
end

function [ss,retcode]=steady_state_evaluation(param_obj,unique_ss,options,ssfunc) %#ok<INUSL>

if nargin<3
    options=optimset('display','none','maxiter',400,'tolfun',1e-6); %#ok<NASGU>
    if nargin<2
        unique_ss=true;
    end
end

retcode=0;
par_vals=vertcat(param_obj.startval);
number_of_regimes=size(par_vals,2);
if unique_ss
    par_vals=mean(par_vals,2);
end
ss=[];
for ii=1:number_of_regimes
    if ii==1 || ~unique_ss
        param=par_vals(:,ii); %#ok<NASGU>
        eval(ssfunc)
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


function lhs=steady_state_residuals(y,func_ss,x,param)
% if the structural model is defined in terms of steady states, then the
% the word ss will appear in the static model. But then, those values are
% the same as those found in y
ss=y;
lhs=func_ss(y,x,ss,param);
end

function lhs=balanced_growth_path_residuals(y,func_bgp,x,param)
ss=y;
% in this case, the y vector also contains the balanced-growth values.
% copying the whole vector to ss does not harm since only the relevant
% elements in the upper part of y will be picked up.
lhs=func_bgp(y,x,ss,param);
end

