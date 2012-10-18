function [T,R,gsig,retcode,optim_opt,itercode,algo_name] = msre_solve(Aminus,A0,Aplus,B,...
    C,sparam,unique_steady_state,...
    Q,order,T0,optim_opt)

% This procedure assumes the steady state has been solved and that apart
% from the constant term representing risk, all variables are all expressed
% as deviations from the steady state

default_solve=struct('solve_fwz_hypothesis',false,...
    'solve_risk',true);
% gather the defaults from fix point iterator as well
default_solve=mergestructures(default_solve,fix_point_iterator());
if nargin==0
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    T=default_solve;
    T.solve_initialization='backward'; %['zeros','zero',{'backward' or 'default'},'random']
    return
end


if nargin<11
    optim_opt=[];
    if nargin<10
        T0=[];
        if nargin<9
            order=1;
            if nargin<8
                error([mfilename,':: at least arguments Aminus,A0,Aplus,B,Q should be passed'])
            end
        end
    end
elseif nargin>11
    error([mfilename,':: number of arguments cannot exceed 8'])
end

%==============================

% just check that the fields are included if not then include them. don't
% go into the business of deleting fields that are irrelevant.
default_fields=fieldnames(default_solve);
for id=1:numel(default_fields)
    if ~isfield(optim_opt,default_fields{id})
        optim_opt.(default_fields{id})=default_solve.(default_fields{id});
    end
end

%==============================
[nn,xx,h]=size(B);

impose_solution=false;
if isempty(T0)
    model_class=any(any(any(Aminus)))+2*any(any(any(Aplus)));
    T0=zeros(size(A0));
    switch model_class
        case {0,2} % static model or forward-looking model
            impose_solution=true;
        case 1 % backward-looking model
            for ii=1:size(A0,3)
                T0(:,:,ii)=-A0(:,:,ii)\Aminus(:,:,ii);
            end
            impose_solution=true;
        case 3 % hybrid model
            %             if strcmp(optim_opt.solve_initialization,'random')
            %                 [junk,junk,junk,sampling_function] = msre_solvers.fwz_newton_system(Aplus,A0,Aminus,Q);
            %             end
            T0=msre_initial_guess(optim_opt.solve_initialization);
    end
end

% if isfield(optim_opt,'accelerate_solver')
%     accelerate=optim_opt.accelerate_solver;
% else
%     accelerate=true;
% end
% acceleration with QR decomposition temporarily removed due to the
% complication arising from averaging the matrices of lead terms in markov
% switching.

% optim_opt=rmfield(optim_opt,'accelerate_solver');

[T,R,gsig,retcode,itercode,algo_name]=...
    markov_switching_rational_expectations_solver(Aminus,A0,Aplus,B,C,T0);

    function T0=msre_initial_guess(solve_initialization)
        T0=zeros(nn,nn,h);
        if h==1
            solve_initialization='zero_init';
            % don't need initialization in the constant-parameter case
        end
        bkwl=any(Aminus(:,:,1),1);
        for istate=2:h
            bkwl=bkwl|any(Aminus(:,:,istate),1);
        end
        switch solve_initialization
            case {'zero','zeros','zero_init'}
            case {'back_init','backward','default'}
                for i_state=1:h
                    % the solution that corresponds to the backward-looking model [default]
                    % this should accelerate the solving and give more accurate results since
                    % we strictly don't need to solve for the non-state columns i.e. the columns
                    % where Aminus is zero
                    T0(:,:,i_state)=-A0(:,:,i_state)\Aminus(:,:,i_state);
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
                            QAT=QAT+A0(:,:,istate);
                            T0(:,:,istate)=-QAT\Aminus(:,:,istate);
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
                n_a0=max(n_a0,norm(A0(:,:,ii_state)));
                n_aplus=max(n_aplus,norm(Aplus(:,:,ii_state)));
                n_aminus=max(n_aminus,norm(Aminus(:,:,ii_state)));
            end
            sn=(n_a0+sqrt(n_a0^2+4*n_aplus*n_aminus))/(2*n_aplus);
       end
    end
    function [T,R,gsig,retcode,itercode,algo_name]=...
            markov_switching_rational_expectations_solver(Aminus,A0,Aplus,B,C,T0)
        % the structural system takes the form
        % E{A_lead(s_{t+1})*X_{t+1}}+A_0(s_t)*X_t+A_lag*X_{t-1}+B*e+C=0.
        % N.B. C is not a constant here. It is an array of first-order derivatives
        % of the system with respect to switching parameters
        % sparam is the matrix of switching parameter values        
                
        if ~isfield(optim_opt,'solver')	|| isempty(optim_opt.solver)
            % normally this will never happen as there is always a default. But if the user inadvertently
            % remove the default, then we have to intervene
            optim_opt.solver=1;
        end
        isevs=is_eigenvalue_solver_candidate(Aplus,A0,Aminus,Q);
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
        
        solve_fwz_hypothesis=optim_opt.solve_fwz_hypothesis;
        solver=optim_opt.solver;
		fwz_flag=false;
        switch solver
            case {0,'msre_klein'}
                iterate_func = @(x)msre_klein(x,Aplus,A0,Aminus,Q,nn,h);
            case {1,'msre_gensys'}
                iterate_func = @(x)msre_gensys(x,Aplus,A0,Aminus,Q,nn,h);
            case {2,'msre_aim'}
                iterate_func = @(x)msre_aim(x,Aplus,A0,Aminus,Q,nn,h);
            case {3,'functional_iteration'}
                iterate_func = @(x)msre_solvers.functional_iteration_(x,Aplus,A0,Aminus,Q,nn,h,solve_fwz_hypothesis);
            case {33,'functional_iteration_'}
                iterate_func = @(x)msre_solvers.functional_iteration(x,Aplus,A0,Aminus,Q,nn,h,solve_fwz_hypothesis);
            case {4,'newton_kronecker'}
                iterate_func = @(x)msre_solvers.newton_kronecker(x,Aplus,A0,Aminus,Q,nn,h,solve_fwz_hypothesis);
            case {44,'newton_kronecker_'}
                iterate_func = @(x)msre_solvers.newton_kronecker_(x,Aplus,A0,Aminus,Q,nn,h,solve_fwz_hypothesis);
            case {5,'newton_system'}
                iterate_func = @(x)msre_solvers.newton_system(x,Aplus,A0,Aminus,Q,nn,h,solve_fwz_hypothesis);
            case {55,'newton_system_'}
                iterate_func = @(x)msre_solvers.newton_system_(x,Aplus,A0,Aminus,Q,nn,h,solve_fwz_hypothesis);
            case {6,'newton_kronecker_iteration'}
                iterate_func = @(x)msre_solvers.newton_kronecker_iteration(x,Aplus,A0,Aminus,Q,nn,h,solve_fwz_hypothesis);
            case {66,'newton_kronecker_iteration_'}
                iterate_func = @(x)msre_solvers.newton_kronecker_iteration_(x,Aplus,A0,Aminus,Q,nn,h,solve_fwz_hypothesis);
			case {7,'fwz_newton_system'}
                [iterate_func,solution_func,inverse_solution_func] = msre_solvers.fwz_newton_system(Aplus,A0,Aminus,Q);
                fwz_flag=true;
            otherwise	% @(x)msre_solvers.(mysolvers{isol})(x,Aplus,A0,Aminus,Q,nn,h,solve_fwz_hypothesis);
                if isnumeric(solver)
                    solver=num2str(solver);
                    error([mfilename,':: solver option ',solver,' unknown'])
                elseif ischar(solver)
                    % if the user has it own solver
                    iterate_func = str2func(solver);
					iterate_func = @(x)iterate_func(x,Aplus,A0,Aminus,Q);
                elseif  isa(solver,'function_handle')
                    iterate_func = @(x)solver(x,Aplus,A0,Aminus,Q);
                else
                    error([mfilename,':: unknown solver option'])
                end
        end
        algo_name=func2str(iterate_func);
        
        % solving for T
        
        solve_risk=optim_opt.solve_risk;
        if impose_solution
            T=T0;
            retcode=0;
            itercode=0;
        else
            if isevs
                [T,junk,retcode]=iterate_func(T0);
                itercode=0;
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
        
        % solving for R and the risk constant
        if retcode
            R=[];
            gsig=[];
        else
            %=============== THIS IS MERELY TO ACCOMMODATE THE OLD MODEL CLASS
            if isempty(unique_steady_state)||size(C,2)<1
                solve_risk=false;
            end
            %=====================
            ATAi=nan(size(T));
            gsig=zeros(nn,h);
            if h>1 && solve_risk
                KK=zeros(nn*h);
                L=zeros(nn,h);
            end
            targ=mean(sparam,2);
            for st=1:h
                st_stretch=(st-1)*nn+1:st*nn;
                QAT=0;
                if h>1 && solve_risk && ~unique_steady_state
                    targ=sparam(:,st);
                end
                %=====================
                for slead=1:h
                    slead_stretch=(slead-1)*nn+1:slead*nn;
                    if solve_fwz_hypothesis
                        QA=Q(st,slead)*Aplus(:,:,st);
                    else
                        QA=Q(st,slead)*Aplus(:,:,slead);
                    end
                    if h>1 && solve_risk
                        KK(st_stretch,slead_stretch)=QA;
                        L(:,st)=L(:,st)+sum(Q(st,slead)*bsxfun(@times,C(:,:,slead),transpose(sparam(:,slead)-targ)),2);
                    end
                    QAT=QAT+QA*T(:,:,slead);
                end
                %=====================
                QAT=QAT+A0(:,:,st);
                if h>1 && solve_risk
                    KK(st_stretch,st_stretch)=KK(st_stretch,st_stretch)+QAT;
                end
                Tst=-QAT\eye(nn);
                ATAi(:,:,st)=Tst;
            end
            if h>1 && solve_risk
                gsig=reshape(-KK\L(:),nn,h);
            end
            
            R=nan(nn,xx,order,h);
            for io=1:order
                for st=1:h
                    if io==1
                        R(:,:,io,st)=ATAi(:,:,st)*B(:,:,st);
                    else
                        Rplus=0;
                        for slead=1:h
                            Rplus=Rplus+Q(st,slead)*R(:,:,io-1,slead);
                        end
                        R(:,:,io,st)=ATAi(:,:,st)*Aplus(:,:,st)*Rplus;
                    end
                end
            end
        end
    end
end

