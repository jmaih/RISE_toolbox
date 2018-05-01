function [sims,retcode]=simulate_nonlinear(obj,varargin)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

% To Do: Conditional forecasting
%---------------------------------
% idea is to condition on future values...
% the system potentially becomes over-identified...
% more conditions than the number of equations to solve for

if isempty(obj)
    
    mydefaults=the_defaults();
    
    if nargout
        
        sims=mydefaults;
        
    else
        
        disp_defaults(mydefaults);
        
    end
    
    return
    
end

sims=[];

obj.options.simul_burn=0;

[obj,retcode]=solve(obj,varargin{:});

if retcode
    
    return
    %     error('model cannot be solved')
end

do_full=~obj.options.simul_sparse;

% initial conditions
%-------------------
Initcond=set_simulation_initial_conditions(obj);

nshocks_periods=size(Initcond.y.econd.data,2);

nshocks=obj.exogenous.number(1);

Initcond.shocks=zeros(nshocks,nshocks_periods); 

Initcond.shocks(Initcond.y.econd.pos,:)=Initcond.y.econd.data(:,:,1);

bigt=obj.options.simul_periods;

% check we have enough periods for the shocks
%---------------------------------------------
missing=bigt-nshocks_periods;

if missing>0
    
    Initcond.shocks=[Initcond.shocks,zeros(nshocks,missing)];
    
end

debug=obj.options.debug;

simul_stack_solve_algo=obj.options.simul_stack_solve_algo;

forecast_exogenous=upload_exogenous_forecast_func();

% without any argument, the function should return the period of the last
% non-zero shock. This is to prevents the algorithms from running
% unnecessary optimizations when the process has already landed.
t_last_shock=forecast_exogenous();

endo_nbr=obj.endogenous.number;

s0=1;

s1=1;

param=obj.parameter_values(:,s1);

ss=obj.solution.ss{s1};

% initial and final conditions
%-----------------------------
sims=ss(:,ones(1,bigt+2)); % initial,simul,final

ystart=Initcond.y.y(obj.inv_order_var,1);%

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

iov=obj.inv_order_var;

ov=obj.order_var;

is_symbolic=strcmp(obj.options.solve_derivatives_type,'symbolic');

is_linear=obj.options.solve_linear;

if is_linear
    
    old_Jac_i=[];
    
    old_Jac=[];
    
end

% sparsity of matrices
%-----------------------
nzlag=nnz(obj.occurrence(:,:,1));

nz0=nnz(obj.occurrence(:,:,2));

nzlead=nnz(obj.occurrence(:,:,3));

nzall=nzlag+nz0+nzlead;

% various indices
%-----------------
[ap_a0_am,np_n0_nm]=symb_jac_reordering_indices(obj.lead_lag_incidence.before_solve,ov,iov);

[numjac_aplus_a0_aminus,a0_ap_1,am_a0_end,am_a0_ap]=...
    num_jac_ordering_indices(obj.lead_lag_incidence.before_solve,do_full);

% locations of various variables in the jacobian
%-----------------------------------------------
[aplus,aplus_cols,a0,aminus,aminus_cols]=get_jacobian_locations();

% location of variables in the system matrices ordered alphabetically
%---------------------------------------------------------------------
[ov_aplus_cols,ov_a0_cols,ov_aminus_cols]=get_alphabetic_system_locations();

[the_leads,the_lags]=dsge_tools.create_endogenous_variables_indices(...
    obj.lead_lag_incidence.before_solve);

the_current=(1:endo_nbr)';

x=[];

if obj.options.simul_recursive
    
    nsmpl=bigt;
    
else
    
    nsmpl=1;
    
end

options=obj.options.optimset;

if debug||obj.options.fix_point_verbose
    
    options.Display='iter';
    
end

options.Jacobian='on';

% options.Algorithm={'levenberg-marquardt',2*.005};

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
                
                [sims(:,t:end),retcode]=do_fsolve(sims(:,t:end));
                
            case 'lsqnonlin'
                
                [sims(:,t:end),retcode]=do_lsqnonlin(sims(:,t:end));

            case 'sparse'
                
                [sims(:,t:end),retcode]=sparse_simulation_algorithm(sims(:,t:end));
                
            case 'lbj'
                
                error('lbj currently broken: post issue on the forum')
                
                [sims(:,t:end),retcode]=relaxation(sims(:,t:end));
                
            otherwise
                % user-defined algorithm
                
                if iscell(simul_stack_solve_algo)
                    
                    user_algo=simul_stack_solve_algo{1};
                    
                    vargs=simul_stack_solve_algo(2:end);
                    
                else
                    
                    vargs={};
                    
                    user_algo=simul_stack_solve_algo;
                    
                end
                
                if ischar(user_algo)
                    
                    user_algo=str2func(user_algo);
                    
                end
                
                x00=vec(sims(:,t+1:end-1));
                
                [x11,retcode]=user_algo(@objective,x00,vargs{:},sims(:,t),...
                    sims(:,end));
                
                sims(:,t+1:end-1)=reshape(x11,endo_nbr,[]);
                                
        end
        
    end
    
end

% store the simulations in a database: use the date for last observation in
% history and not first date of forecast
%--------------------------------------------------------------------------

if isempty(obj.options.simul_history_end_date)
    
    obj.options.simul_history_end_date=0;
    
end

start_date=date2serial(Initcond.simul_history_end_date);%-y0cols+1

sims=ts(start_date,full(sims)',obj.endogenous.name);

sims=pages2struct(sims);

    function [sims,retcode]=do_lsqnonlin(sims)
        
        sims=reshape(sims,endo_nbr,[]);
        
        Y0=sims(:,2:end-1);
        
        [y1,~,~,exitflag]=lsqnonlin(@objective,Y0(:),[],[],options,...
            sims(:,1),sims(:,end));
        
        if exitflag==0
            
            retcode=21;
            
        elseif ~utils.error.valid(y1)
            
            retcode=22;
            
        else
            
            retcode=(~ismember(exitflag,[1,2,3,4]))*23;
            
        end
        
        if retcode==0
            
            sims(:,2:end-1)=reshape(y1,endo_nbr,[]);
            
        end
        
    end

    function [sims,retcode]=do_fsolve(sims)
        
        Y0=sims(:,2:end-1);
        
        [y1,~,exitflag]=fsolve(@objective,Y0(:),options,sims(:,1),...
            sims(:,end));
        
        if exitflag==0
            
            retcode=21;
            
        elseif ~utils.error.valid(y1)
            
            retcode=22;
            
        else
            
            retcode=(~ismember(exitflag,[1,2,3,4]))*23;
            
        end
        
        if retcode==0
            
            sims(:,2:end-1)=reshape(y1,endo_nbr,[]);
            
        end
        
    end

    function [BIGX,retcode]=relaxation(BIGX)
        
        BIGX=reshape(BIGX,endo_nbr,[]);
        
        nsyst=size(BIGX,2)-2;
        
        [BIGX,~,retcode]=fix_point_iterator(@relaxation_iterator,BIGX,...
            obj.options);
        
        function [BIGX,crit]=relaxation_iterator(BIGX)
            
            Y=build_lead_current_lag_matrix(BIGX);
            
            d=residuals(dynamic_model,Y,...
                x(:,1:nsyst),ss,param,def,s0,s1);
            
            C=cell(nsyst,1);
            
            B=C;
            
            A=C;
            
            D=C;
            
            relaxed_jacobian();
            
            for isyst=1:nsyst
                
                if isyst==1
                    
                    D{isyst}=B{isyst}\C{isyst};
                    
                    d(:,isyst)=B{isyst}\d(:,isyst);
                    
                else
                    
                    tmp=(B{isyst}-A{isyst}*D{isyst-1})\eye(endo_nbr);
                    
                    if isyst<nsyst
                        
                        D{isyst}=tmp*C{isyst};
                        
                    end
                    
                    d(:,isyst)=tmp*(d(:,isyst)+A{isyst}*d(:,isyst-1));
                    
                end
                
            end
            
            dx=BIGX(:,2:end-1);
            
            for isyst=nsyst:-1:1
                
                if isyst==nsyst
                    
                    dx(:,isyst)=-d(:,isyst);
                    
                else
                    
                    dx(:,isyst)=-d(:,isyst)-D{isyst}*d(:,isyst+1);
                    
                end
                
            end
            
            crit=full(max(max(abs(dx))));
            
            if nsyst==1
                
                dx=dx(:,1);
                
            end
            
            s=junior_maih_2016(BIGX,dx,nsyst);
            
            BIGX(:,2:end-1)=BIGX(:,2:end-1)+s*dx;
            
            function relaxed_jacobian()
                
                for isyst_=1:nsyst
                    
                    [C{isyst_},B{isyst_},A{isyst_}]=...
                        build_system_matrices(Y(:,isyst_),isyst_);
                    
                end
                
            end
            
        end
        
    end

    function [BIGX,retcode]=sparse_simulation_algorithm(BIGX)
        
        BIGX=reshape(BIGX,endo_nbr,[]);
        
        nsyst=size(BIGX,2)-2;
        
        [BIGX,~,retcode]=fix_point_iterator(@sparse_iterator,BIGX,...
            obj.options);
        
        function [BIGX,crit]=sparse_iterator(BIGX)
                        
            Y=build_lead_current_lag_matrix(BIGX);
            
            resid=residuals(dynamic_model,Y,...
                x(:,1:nsyst),ss,param,def,s0,s1);
            
            if is_linear
                
                if isempty(old_Jac_i)
                    
                    A=build_big_jacobian(Y,nsyst);
                    
                    old_Jac_i=A\eye(size(A,1));
                    
                end
                
                dx=reshape(-old_Jac_i*resid(:),endo_nbr,[]);
                
            else
                
                A=build_big_jacobian(Y,nsyst);
                
                dx=reshape(-A\resid(:),endo_nbr,[]);
                
            end
            
            crit=full(max(max(abs(dx))));
            
            if nsyst==1
                
                dx=dx(:,1);
                
            end
                        
            s=junior_maih_2016(BIGX,dx,nsyst);
            
            BIGX(:,2:end-1)=BIGX(:,2:end-1)+s*dx;
            
        end
        
    end

    function s=junior_maih_2016(BIGX,dx,nsyst)
        % Now we need to ensure that the step is valid
        
        s=1;
        
        if ~utils.error.valid(dx)
                        
            return
            
        end
        
        iter=0;
        
        while s>0 && ~utils.error.valid(update_residuals(s))
            
            s=0.9*s;
            
            iter=iter+1;
            
            if debug
                
                fprintf(1,'iteration %0.0f, slowc %0.8f\n',iter,s);
                
            end
            
        end
        
        if s==0
            % set s to nan in order to avoid an "infinite" (maxiter) loop
            
            s=nan;
            
        end
        
        function resid1=update_residuals(s)
            
            BIGX1=BIGX;
            
            BIGX1(:,2:end-1)=BIGX1(:,2:end-1)+s*dx;
            
            Y1=build_lead_current_lag_matrix(BIGX1);
            
            resid1=residuals(dynamic_model,Y1,...
                x(:,1:nsyst),ss,param,def,s0,s1);
            
        end
        
    end

    function [ov_aplus_cols,ov_a0_cols,ov_aminus_cols]=get_alphabetic_system_locations()
                
        ov_aplus_cols=ov(aplus_cols);
        
        ov_aminus_cols=ov(aminus_cols);
        
        ov_a0_cols=ov;

    end
            
    function [aplus,aplus_cols,a0,aminus,aminus_cols]=get_jacobian_locations()
        
        ns=numel(obj.locations.before_solve.v.s_0);
        
        np=numel(obj.locations.before_solve.v.p_0);
        
        nb=numel(obj.locations.before_solve.v.b_0);
        
        nf=numel(obj.locations.before_solve.v.f_0);
        
        aplus_cols=ns+np+(1:nb+nf);
        
        aplus=obj.locations.before_solve.v.bf_plus;
        
        aminus_cols=ns+(1:np+nb);
        
        aminus=obj.locations.before_solve.v.pb_minus;
        
        a0=obj.locations.before_solve.v.t_0;
        
    end

    function [resid,Jac]=objective(X,X0,Xlast)
        
        X=reshape(X,endo_nbr,[]);
        
        nsyst=size(X,2);
        
        X=[X0,X,Xlast];
        
        Y=build_lead_current_lag_matrix(X);
        
        resid=residuals(dynamic_model,Y,x(:,1:nsyst),...
            ss,param,def,s0,s1);
        
        resid=resid(:);
        
        if nargout>1
            
            if is_linear
                
                if isempty(old_Jac)
                    
                    old_Jac=build_big_jacobian(Y,nsyst);
                    
                end
                
                Jac=old_Jac;
                
            else
                
                Jac=build_big_jacobian(Y,nsyst);
                
            end
            
        end
        
    end

    function J=build_big_jacobian(Y,nsyst)
        
        if nargin<2
            
            nsyst=size(Y,2);
            
        end
        
        nrows=nsyst*endo_nbr;
        
        ncols=nrows;
        
        nzmax=max_non_zeros(nsyst);
        
        J=spalloc(nrows,ncols,nzmax);
        % J=sparse([],[],[],nrows,ncols,nzmax);
        
        offset_rows=0;
        
        offset_cols=0;
        
        for isyst=1:nsyst
            
            [Aplus,A0,Aminus]=build_system_matrices(Y(:,isyst),isyst);
            
            if isyst>2
                
                offset_cols=offset_cols+endo_nbr;
                
            end
                
            if isyst==nsyst
                
                J(offset_rows+(1:endo_nbr),offset_cols+am_a0_end)=[Aminus,A0];
                        
            elseif isyst==1
                
                J(1:endo_nbr,offset_cols+a0_ap_1)=[A0,Aplus];
                
            else
                
                J(offset_rows+(1:endo_nbr),offset_cols+am_a0_ap)=[Aminus,A0,Aplus];
                
            end
            
            offset_rows=offset_rows+endo_nbr;
            
        end
        
    end

    function n=max_non_zeros(nsyst)
        
        n=nzall*(nsyst-2)+(nz0+nzlead)+(nzlag+nz0);
        
    end

    function Y=build_lead_current_lag_matrix(X)
        
        Y=[
            X(the_leads,3:end)
            X(the_current,2:end-1)
            X(the_lags,1:end-2)
            ];
        
    end

    function [Aplus,A0,Aminus]=build_system_matrices(y,isyst)
        
        if is_symbolic
            
            SymbJac=utils.code.evaluate_functions(jacobian_dynamic_model,...
                y,x(:,isyst),ss,param,def,s0,s1);
            
            [Aplus,A0,Aminus]=decompose_symbolic_jacobian(SymbJac);
            
        else
            
            NumJac=numerical_jacobian(dynamic_model,y,x(:,isyst),...
                ss,param,def,s0,s1);
            
            [Aplus,A0,Aminus]=decompose_numerical_jacobian(NumJac);
            
        end
        
        if debug
            if is_symbolic
                % compute numeric
                
                NumJac=numerical_jacobian(dynamic_model,y,x(:,isyst),...
                    ss,param,def,s0,s1);
                
                [Aplus_,A0_,Aminus_]=decompose_numerical_jacobian(NumJac);
            else
                % compute symbolic
                SymbJac=utils.code.evaluate_functions(jacobian_dynamic_model,...
                    y,x(:,isyst),ss,param,def,s0,s1);
                
                [Aplus_,A0_,Aminus_]=decompose_symbolic_jacobian(SymbJac);
            end
            
            max(max(abs(Aplus_-Aplus)))
            
            max(max(abs(A0_-A0)))
            
            max(max(abs(Aminus_-Aminus)))
            
        end
        
    end

    function [Aplus,A0,Aminus]=decompose_symbolic_jacobian(Jac)
        
        if do_full
            
            Aplus=zeros(endo_nbr);
            
            Aminus=Aplus;
            
            A0=Aplus;
            
            Aplus(:,ov_aplus_cols)=Jac(:,aplus);
            
            Aminus(:,ov_aminus_cols)=Jac(:,aminus);
            
            A0(:,ov_a0_cols)=Jac(:,a0);
            
        else
            
            Aplus=Jac(:,ap_a0_am{1});
            
            A0=Jac(:,ap_a0_am{2});
            
            Aminus=Jac(:,ap_a0_am{3});
            
        end
        
    end

    function [Aplus,A0,Aminus]=decompose_numerical_jacobian(Jac)
        
        if do_full
            
            biga=zeros(endo_nbr,3*endo_nbr);
            
            biga(:,numjac_aplus_a0_aminus)=Jac;
            
            Aplus=biga(:,1:endo_nbr);
            
            A0=biga(:,endo_nbr+1:2*endo_nbr);
            
            Aminus=biga(:,2*endo_nbr+1:end);
            
        else
            
            aplus_=1:np_n0_nm(1);
            
            a0_=np_n0_nm(1)+(1:np_n0_nm(2));
            
            aminus_=np_n0_nm(1)+np_n0_nm(2)+(1:np_n0_nm(3));
            
            Aplus=Jac(:,aplus_);
            
            A0=Jac(:,a0_);
            
            Aminus=Jac(:,aminus_);
            
        end
        
    end

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

    function resid=residuals(varargin)
        % due to a change in behavior in evaluate_functions, a special
        % reshaping is needed here. If the behavior is undone, changing
        % this will be necessary as well.
        resid=utils.code.evaluate_functions(varargin{:});
        
% % % % %         resid=reshape(resid,[],endo_nbr).';
        
    end

end

function [numjac_aplus_a0_aminus,a0_ap_1,am_a0_end,am_a0_ap]=...
    num_jac_ordering_indices(lli,do_full)

n=size(lli,1);

aplus=find(lli(:,1));

aplus=aplus(:).';

np=numel(aplus);

a0=(1:n);

aminus=find(lli(:,3));

aminus=aminus(:).';

numjac_aplus_a0_aminus=[aplus,np+a0,np+n+aminus];

if do_full
    
    a0_ap_1=1:2*n;% =[A0,Aplus];
    
    am_a0_end=1:2*n;%=[Aminus,A0];
    
    am_a0_ap=1:3*n;%=[Aminus,A0,Aplus];
               
else
    
    a0_ap_1=[a0,n+aplus];% =[A0,Aplus];
    
    am_a0_end=[aminus,n+a0];%=[Aminus,A0];
    
    am_a0_ap=[am_a0_end,2*n+aplus];%=[Aminus,A0,Aplus];
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

function [ap_a0_am,np_n0_nm]=symb_jac_reordering_indices(lli,ov,iov)

offset=0;

np_n0_nm=sum(lli>0,1);

ap_a0_am=cell(1,3);

for icol=1:3
    
    if icol==2
        % current
        %---------
        
        batch=np_n0_nm(icol-1)+iov;
        
    else
        % leads and lags
        %----------------
        
        l=lli(:,icol);
        
        lov=l(ov);
        
        lf=l(l>0);
        
        lovf=lov(lov>0);
        
        batch=zeros(1,np_n0_nm(icol));
        
        for jj=1:np_n0_nm(icol)
            
            tmp=find(lf(jj)==lovf);
            
            batch(jj)=offset+tmp;
            
        end
        
    end
    
    offset=offset+np_n0_nm(icol);
    
    ap_a0_am{icol}=batch;
    
end

end

function d=the_defaults()

d={'simul_exogenous_func','',@(x)isa(x,'function_handle'),...
    'simul_exogenous_func must be a function handle' 
    
    'simul_stack_solve_algo','sparse',...
    @(x)ismember(x,{'fsolve','lsqnonlin','sparse','lbj'}),...
    'simul_stack_solve_algo must be ''fsolve'',''lsqnonlin'',''sparse'' or ''lbj'''
    
    'simul_recursive',false,@(x)islogical(x),'simul_recursive must be a logical'
    
    'simul_sparse',true,@(x)islogical(x),'simul_sparse must be a logical'
    };
 
end