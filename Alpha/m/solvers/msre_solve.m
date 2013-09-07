function [T,R,gsig,theta_hat,retcode,optim_opt,itercode,algo_name] = msre_solve(Aminus,A0,Gplus01,B,...
    Gtheta,sparam,unique_steady_state,...
    Q,T0,optim_opt)

% This procedure assumes the steady state has been solved and that apart
% from the constant term representing risk, all variables are all expressed
% as deviations from the steady state
solve_disable_theta=true;

if nargin==0
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    default_solve=struct('solve_risk',true,...
        'solve_initialization','backward',...
        'solve_disable_theta',true); % 'solve_fwz_hypothesis',false,...
    % gather the defaults from fix point iterator as well
    default_solve=mergestructures(default_solve,fix_point_iterator());
    T=default_solve;
    %solve_initialization=['zeros','zero',{'backward' or 'default'},'random']
    return
end

if nargin~=10
    error([mfilename,':: number of arguments should be 10'])
end
solve_disable_theta=optim_opt.solve_disable_theta;
%==============================
% just check that the fields are included if not then include them. don't
% go into the business of deleting fields that are irrelevant.
solve_expect_order=optim_opt.solve_expect_order;
%==============================
h=numel(B);
nn=size(B{1},1);

impose_solution=false;
if isempty(T0)
    M=false;
    P=false;
    for istate0=1:h
        M=M||any(Aminus{istate0}(:));
        for istate1=1:h
            P=P||any(Gplus01{istate0,istate1}(:));
        end
    end
    model_class=M+2*P;
%     warnstate=warning('query','all');
%     warning('off','MATLAB:singularMatrix')
    switch model_class
        case {0,2} % static model or forward-looking model
            T0=msre_initial_guess('zeros');
            impose_solution=true;
        case 1 % backward-looking model
            T0=msre_initial_guess('backward');
            impose_solution=true;
        case 3 % hybrid model
            if h==1 && model_class==3
                T0=msre_initial_guess('zero_init');
            else
                T0=msre_initial_guess(optim_opt.solve_initialization);
            end
    end
%     warning(warnstate)
end

[T,R,gsig,theta_hat,retcode,itercode,algo_name]=...
    markov_switching_rational_expectations_solver(Aminus,A0,Gplus01,B,Gtheta,T0);

    function T0=msre_initial_guess(solve_initialization)
        T0=zeros(nn,nn,h);
        bkwl=any(Aminus{1},1);
        for istate=2:h
            bkwl=bkwl|any(Aminus{istate},1);
        end
        switch solve_initialization
            case {'zero','zeros','zero_init'}
            case {'back_init','backward','default'}
                for i_state=1:h
                    % the solution that corresponds to the backward-looking model [default]
                    % this should accelerate the solving and give more accurate results since
                    % we strictly don't need to solve for the non-state columns i.e. the columns
                    % where Aminus is zero
                    T0(:,:,i_state)=-A0{i_state}\Aminus{i_state};
                    if any(any(isnan(T0(:,:,i_state))))
                        T0(:,:,i_state)=-pinv(A0{i_state})*Aminus{i_state};
                        % the A0 matrix can be singular, for instance in the
                        % case of the zlb
                    end
                end
            case {'random','rand_init'}
                level=3;
                switch level
                    case 0
                        T0(:,bkwl,:)=randn(nn,sum(bkwl),h);
                    case 1
                        AT=zeros(nn,nn,h);
                        sn=solvent_norm();
                        AT(:,bkwl,:)=sn^2*randn(nn,sum(bkwl),h);
                        for istate=1:h
                            QAT=zeros(nn);
                            for jstate=1:h
                                QAT(:,bkwl)=QAT(:,bkwl)+Q(istate,jstate)*AT(:,bkwl,jstate);
                            end
                            QAT=QAT+A0{istate};
                            T0(:,:,istate)=-QAT\Aminus{istate};
                        end
                    case 2
                        Tbkw=msre_initial_guess('back_init');
                        T0(:,bkwl,:)=Tbkw(:,bkwl,:).*rand(nn,sum(bkwl),h);
                    case 3
                        sn=solvent_norm();
                        T0(:,bkwl,:)=sn*randn(nn,sum(bkwl),h);
                end
        end
        function sn=solvent_norm()
            n_a0=0;n_aplus=0;n_aminus=0;
            for ii_state=1:h
                n_a0=max(n_a0,norm(A0{ii_state}));
                aplus=0;
                for s1=1:h
                    aplus=aplus+Gplus01{ii_state,s1};
                end
                n_aplus=max(n_aplus,norm(aplus));
                n_aminus=max(n_aminus,norm(Aminus{ii_state}));
            end
            sn=(n_a0+sqrt(n_a0^2+4*n_aplus*n_aminus))/(2*n_aplus);
        end
    end
    function [T,R,msig,theta_hat,retcode,itercode,algo_name]=...
            markov_switching_rational_expectations_solver(Aminus,A0,Gplus01,B,Gtheta,T0)
        % the structural system takes the form
        % E{A_lead(s_{t+1})*X_{t+1}}+A_0(s_t)*X_t+A_lag*X_{t-1}+B*e+Gtheta=0.
        % N.B. Gtheta is not a constant here. It is an array of first-order derivatives
        % of the system with respect to switching parameters
        % sparam is the matrix of switching parameter values
        
        % this to ensure that the constant-parameter case does not crash,
        % we initialize theta_hat in the following way
        theta_hat=repmat({zeros(0,1)},1,h); % <-- cell(1,h) 
        
        if ~isfield(optim_opt,'solver')	|| isempty(optim_opt.solver)
            % normally this will never happen as there is always a default. But if the user inadvertently
            % remove the default, then we have to intervene
            optim_opt.solver=1;
        end
        isevs=is_eigenvalue_solver_candidate(Gplus01,A0,Aminus,Q);
        if ischar(optim_opt.solver)
            if ismember(optim_opt.solver,{'msre_gensys','msre_aim','msre_klein','functional_iteration',...
                    'newton_kronecker','newton_system','newton_kronecker_iteration'})
                if isevs && ~ismember(optim_opt.solver,{'msre_gensys','msre_aim','msre_klein'})
                    optim_opt.solver=0;	% Klein by default
                elseif ~isevs && ismember(optim_opt.solver,{'msre_gensys','msre_aim','msre_klein'})
                    optim_opt.solver=3; % functional iteration by default
                end
            else
                % then the user is providing a function that is not implemented in the toolbox and that's ok as long as it respects the syntax
            end
        elseif isnumeric(optim_opt.solver)
            if ~ismember(optim_opt.solver,(0:6))
                error([mfilename,':: option ',num2str(optim_opt.solver),' not in the range [0,6]'])
            end
            if isevs && ~ismember(optim_opt.solver,[0,1,2])
                optim_opt.solver=1;	% gensys by default
            elseif ~isevs && ismember(optim_opt.solver,[0,1,2])
                optim_opt.solver=3; % functional iteration by default
            end
        end
        %        solve_fwz_hypothesis=optim_opt.solve_fwz_hypothesis;
        solver=optim_opt.solver;
        fwz_flag=false;
        switch solver
            case {0,'msre_klein'}
                iterate_func = @(x)msre_klein(x,Gplus01,A0,Aminus,Q,nn,h);
            case {1,'msre_gensys'}
                iterate_func = @(x)msre_gensys(x,Gplus01,A0,Aminus,Q,nn,h);
            case {2,'msre_aim'}
                iterate_func = @(x)msre_aim(x,Gplus01,A0,Aminus,Q,nn,h);
            case {3,'functional_iteration'}
                iterate_func = @(x)msre_solvers.functional_iteration(x,Gplus01,A0,Aminus,Q,nn,h);
            case {4,'newton_kronecker'}
                iterate_func = @(x)msre_solvers.newton_kronecker(x,Gplus01,A0,Aminus,Q,nn,h);
            case {5,'newton_system'}
                iterate_func = @(x)msre_solvers.newton_system(x,Gplus01,A0,Aminus,Q,nn,h);
            case {6,'newton_kronecker_iteration'}
                iterate_func = @(x)msre_solvers.newton_kronecker_iteration(x,Gplus01,A0,Aminus,Q,nn,h);
            case {7,'fwz_newton_system'}
                [iterate_func,solution_func,inverse_solution_func] = msre_solvers.fwz_newton_system(Gplus01,A0,Aminus,Q);
                fwz_flag=true;
            otherwise	% @(x)msre_solvers.(mysolvers{isol})(x,Gplus01,A0,Aminus,Q,nn,h);
                if isnumeric(solver)
                    solver=num2str(solver);
                    error([mfilename,':: solver option ',solver,' unknown'])
                elseif ischar(solver)
                    % if the user has it own solver
                    iterate_func = str2func(solver);
                    iterate_func = @(x)iterate_func(x,Gplus01,A0,Aminus,Q);
                elseif  isa(solver,'function_handle')
                    iterate_func = @(x)solver(x,Gplus01,A0,Aminus,Q);
                else
                    error([mfilename,':: unknown solver option'])
                end
        end
        algo_name=func2str(iterate_func);
        
        % solving for T
        
        if impose_solution
            T=T0;
            retcode=0;
            itercode=0;
        else
            if isevs
                [T,~,retcode]=iterate_func(T0);
                itercode=0;
                if ~retcode
                    T{1}=sparse(T{1});
                end
            else
                if fwz_flag % transform T into X
                    T0=inverse_solution_func(T0);
                end
                [T,itercode,retcode]=fix_point_iterator(iterate_func,T0,optim_opt);
                % N.B.: the last option is special for switching solvers
                if fwz_flag && ~retcode % transform X into T
                    T=solution_func(T);
                end
            end
        end
        if ~retcode && ~iscell(T)
            % put in cell form
            tmp=cell(1,h);
            for istate=1:h
                tmp{istate}=sparse(T(:,:,istate));
            end
            T=tmp; clear tmp
        end
        
        if retcode
            R=[];
            msig=[];
        else
            A0Gm_inv=msre_solvers.compute_preconditioner(Gplus01,A0,Q,T);
            R=cell(1,h);
            for io=1:solve_expect_order
                for st=1:h
                    if io==1
                        R{st}(:,:,io)=-A0Gm_inv{st}*B{st};
                    else
                        Rplus=0;
                        for slead=1:h
                            Rplus=Rplus+Q(st,slead)*R{slead}(:,:,io-1);
                        end
                        R{st}(:,:,io)=-A0Gm_inv{st}*Gplus01{st}*Rplus;
                    end
                end
            end
            
            % solving for R and the risk constant
            %------------------------------------
            solve_risk=h>1 && optim_opt.solve_risk && size(Gtheta{1,1},2)>=1;
            msig=zeros(nn,h);
            if solve_risk
                % precondition
                %-------------
                Gbarplus01=Gplus01;
                Dbartheta=zeros(nn,h);
                if unique_steady_state
                    target=0;
                    [PAI00,retcode]=initial_markov_distribution(Q,'ergodic');
                    if retcode
                        return
                    end
                    for s0=1:h
                        target=target+PAI00(s0)*sparam(:,s0);
                    end
                end
                for s0=1:h
                    if ~unique_steady_state
                        target=sparam(:,s0);
                    end
                    for s1=1:h
                        if s0==1
                            theta_hat{s1}=sparam(:,s1)-target;
                        end
                        Gbarplus01{s0,s1}=A0Gm_inv{st}*Gplus01{s0,s1};
                        Dbartheta(:,s0)=Dbartheta(:,s0)+...
                            A0Gm_inv{st}*Gtheta{s0,s1}*theta_hat{s1};
                    end
                end
                clear Gplus01 Gtheta
                   
                msigfunc=@(x)m_sig_matrix_vector(x,Gbarplus01,nn,h);
                [msig,retcode,tau]=transpose_free_quasi_minimum_residual(...
                    msigfunc,... % coefficient matrix
                    -(1-solve_disable_theta)*Dbartheta(:),... % right hand side
                    msig(:),... % initial guess
                    optim_opt.fix_point_TolFun,... % tolerance level
                    optim_opt.fix_point_maxiter,... % maximum number of iterations
                    optim_opt.fix_point_verbose); % flag for printing progress or not
                msig=reshape(msig,nn,h);
                if retcode==201 && tau<=optim_opt.fix_point_TolFun
                    retcode=0;
                end
            end
            if retcode==0
                tmp=msig;
                msig=cell(1,h);
                for st=1:h
                    msig{st}=tmp(:,st);
                end
            end
        end
    end
end

function Ax=m_sig_matrix_vector(msig,Gbarplus01,n,h)
Ax=zeros(n,h);
msig=reshape(msig,n,h);
for s0=1:h
    Ax(:,s0)=msig(:,s0);
    for s1=1:h
        % msig(:,s0)+ sum_s1(Q(s0,s1)*Aplus{s0}*msig(:,s1));
        Ax(:,s0)=Ax(:,s0)+Gbarplus01{s0,s1}*msig(:,s1);
    end
end
Ax=Ax(:);
end

