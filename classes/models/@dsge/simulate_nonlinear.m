function [sims,retcode]=simulate_nonlinear(obj,varargin)
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
    sims=struct('simul_stack_solve_algo','sparse',...
        'simul_recursive',false,...
        'simul_exogenous_func','');
    return
end
sims=[];

obj.options.simul_burn=0;

[obj,retcode]=solve(obj,varargin{:});

if retcode
    return
    %     error('model cannot be solved')
end

% initial conditions
%-------------------
Initcond=set_simulation_initial_conditions(obj);
bigt=obj.options.simul_periods;

simul_stack_solve_algo=obj.options.simul_stack_solve_algo;

forecast_exogenous=upload_exogenous_forecast_func();

% without any argument, the function should return the period of the last
% non-zero shock. This is to prevents the algorithms from running
% unnecessary optimizations when the process has already landed.
t_last_shock=forecast_exogenous();

endo_nbr=obj.endogenous.number(end);

s0=1;
s1=1;
param=obj.parameter_values(:,s1);
sparam=param;
ss=obj.solution.ss{s1};
% initial and final conditions
%-----------------------------
sims=ss(:,ones(1,bigt+2)); % initial,simul,final
ystart=Initcond.y.y;
yfinal=ss;
sims(:,1)=ystart;
sims(:,end)=yfinal;

def=obj.solution.definitions{s1};
% dynamic model
%--------------
dynamic_model=obj.routines.dynamic;
if isa(obj.routines.probs_times_dynamic_derivatives,'function_handle')
    jacobian_dynamic_model=@(varargin)obj.routines.probs_times_dynamic_derivatives(varargin{:});
else
    order__=1;
    jacobian_dynamic_model=obj.routines.probs_times_dynamic_derivatives(order__);
end
order_var=obj.order_var.after_solve;
% locations of various variables in the jacobian
%-----------------------------------------------
ns=numel(obj.locations.before_solve.v.s_0);
np=numel(obj.locations.before_solve.v.p_0);
nb=numel(obj.locations.before_solve.v.b_0);
nf=numel(obj.locations.before_solve.v.f_0);
aplus_cols=ns+np+(1:nb+nf);
aplus=obj.locations.before_solve.v.bf_plus;
aminus_cols=ns+(1:np+nb);
aminus=obj.locations.before_solve.v.pb_minus;
a0=obj.locations.before_solve.v.t_0;
jac_aminus_a0_aplus=[aminus,a0,aplus];
jac_aminus_a0=[aminus,a0];
jac_a0_aplus=[a0,aplus];
A_aminus_a0=[aminus_cols,endo_nbr+(1:endo_nbr)];
A_aminus_a0_aplus=[A_aminus_a0,2*endo_nbr+aplus_cols];
A_a0_aplus=[1:endo_nbr,endo_nbr+aplus_cols];

lead_lag_incidence=obj.lead_lag_incidence.after_solve;
the_leads=find(lead_lag_incidence(:,1)>0);
the_current=(1:endo_nbr)';
the_lags=find(lead_lag_incidence(:,3)>0);
exo_nbr=sum(obj.exogenous.number);
x=[];
if obj.options.simul_recursive
    nsmpl=bigt;
else
    nsmpl=1;
end
for t=1:min(t_last_shock,nsmpl)
    if ~retcode
        % compute the forecast of the exogenous
        %--------------------------------------
        x=forecast_exogenous(obj,t);
        
        % initialize the solver at the steady state
        %-------------------------------------------
        fprintf(1,'\n pass # %0.0f of %0.0f\n',t,nsmpl);
        switch simul_stack_solve_algo
            case 'fsolve'
                options=obj.options.optimset;
                if obj.options.debug||obj.options.fix_point_verbose
                    options.Display='iter';
                end
                Y0=sims(:,1+(t:bigt));
                [sims(:,1+(t:bigt)),~,exitflag]=...
                    fsolve(@objective,Y0,options,sims(:,t),sims(:,end));
                %-------------------------------------------
                if exitflag==0
                    retcode=21;
                elseif ~utils.error.valid(sims)
                    retcode=22;
                else
                    retcode=(~ismember(exitflag,[1,2,3,4]))*23;
                end
                %-------------------------------------------
            case 'sparse'
                if t==100
                    keyboard
                end
                [sims(:,t:end),retcode]=sparse_simulation_algorithm(sims(:,t:end));
            case 'lbj'
                error('LBJ''s relaxation algorithm not yet implemented')
            otherwise
                error(['unknown algorithm ',parser.any2str(simul_stack_solve_algo)])
        end
    end
end
% store the simulations in a database: use the date for last observation in
% history and not first date of forecast
%--------------------------------------------------------------------------
if isempty(obj.options.simul_history_end_date)
    obj.options.simul_history_end_date=0;
end
start_date= obj.options.simul_history_end_date;
sims=ts(start_date,full(sims)',obj.endogenous.name);
sims=pages2struct(sims);


    function forecast_exogenous=upload_exogenous_forecast_func()
        
        forecast_exogenous=obj.options.simul_exogenous_func;
        if isempty(forecast_exogenous)
            forecast_exogenous=@x_forecast;
        end
        function x=x_forecast(~,t)
            if nargin==0
                last_shock_period=find(any(Initcond.shocks,1),1,'last');
                x=last_shock_period;
            else
                if obj.options.simul_recursive
                    x=zeros(size(Initcond.shocks));
                    x(:,1)=Initcond.shocks(:,t);
                else
                    x=Initcond.shocks(:,1:end);
                end
            end
        end
    end

    function resid=objective(X,X0,Xlast)
        X=reshape(X,endo_nbr,[]);
        nsyst=size(X,2);
        resid=X;
        X=[X0,X,Xlast];
        for isyst=1:nsyst
            y=build_lead_current_lag_vector(X,isyst);
            resid(:,isyst)=utils.code.evaluate_functions(dynamic_model,...
                y,x(:,isyst),ss,param,sparam,def,s0,s1);
        end
    end

    function y=build_lead_current_lag_vector(X,iter)
        XX=X(:,iter:iter+2);
        
        y=[ % indices=[the_leads;the_current;the_lags];
            XX(the_leads,end)
            XX(the_current,end-1)
            XX(the_lags,end-2)
            ];
    end

    function [Y,retcode]=sparse_simulation_algorithm(Y)
        Y=reshape(Y,endo_nbr,[]);
        nsyst=size(Y,2)-2;
        A=initialize_jacobian(nsyst);
        resid=nan(endo_nbr,nsyst);
        
        [Y,~,retcode]=fix_point_iterator(@iteration_engine,Y,obj.options);
        
        function [Y,crit]=iteration_engine(Y)
            irows=1:endo_nbr;
            offset_cols=0;
            for isyst=1:nsyst
                y=build_lead_current_lag_vector(Y,isyst);%,nsyst
                
                resid(:,isyst)=utils.code.evaluate_functions(dynamic_model,...
                    y,x(:,isyst),ss,param,sparam,def,s0,s1);
                
                SymbJac=utils.code.evaluate_functions(jacobian_dynamic_model,...
                    y,x(:,isyst),ss,param,sparam,def,s0,s1);
                if obj.options.debug
                    NumJac=numerical_jacobian(dynamic_model,y,x(:,isyst),ss,param,...
                        sparam,def,s0,s1);
                    NumJac0=NumJac(:,nb+nf+(1:endo_nbr));
                    % the jacobian comes in the order_var order. I have to
                    % re-order the NumJac portion to test
                    tmp=NumJac0(:,order_var)-full(SymbJac(:,nb+nf+(1:endo_nbr)));
                    crit0=max(abs(tmp(:)));
                    if crit0>1e-4
                        keyboard
                    end
                    fprintf(1,'%0.8f\n',crit0);
                end
                
                if isyst==1
                    A(irows,A_a0_aplus)=SymbJac(:,jac_a0_aplus);
                elseif isyst==nsyst
                    A(irows,offset_cols+A_aminus_a0)=SymbJac(:,jac_aminus_a0);
                else
                    A(irows,offset_cols+A_aminus_a0_aplus)=SymbJac(:,jac_aminus_a0_aplus);
                end
                
                if isyst<nsyst
                    irows=irows+endo_nbr;
                    if isyst>=2
                        offset_cols=offset_cols+endo_nbr;
                    end
                end
            end
            dx=-A\resid(:);
            
            crit=max(abs(dx));
            dx=reshape(dx,endo_nbr,[]);
            if nsyst==1
                dx=dx(:,1);
            end
            Y(order_var,2:end-1)=Y(order_var,2:end-1)+dx;
        end
    end

    function A=initialize_jacobian(nperiods)
        y=[ % indices=[the_leads;the_current;the_lags];
            ss(the_leads,1)
            ss(the_current,1)
            ss(the_lags,1)
            ];
        ll0_nbr=numel(y);
        JJ=utils.code.evaluate_functions(jacobian_dynamic_model,...
            y,rand(exo_nbr,1),ss,param,sparam,def,s0,s1);
        JJ=JJ(:,1:ll0_nbr);
        nzmax=nnz(JJ)*nperiods;
        nrows=endo_nbr*nperiods;
        ncols=nrows;
        A=sparse([],[],[],nrows,ncols,nzmax);
    end
end

function J=numerical_jacobian(funcs,y,varargin)
nrows=numel(funcs);
ncols=numel(y);
J=nan(nrows,ncols);
for irow=1:nrows
    J(irow,:)=utils.numdiff.jacobian(funcs{irow},y,varargin{:});
end
end