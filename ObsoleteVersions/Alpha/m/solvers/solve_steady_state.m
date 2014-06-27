function [ys,retcode]=solve_steady_state(ys0,x_ss,~,param,def,~,~,...
    resid_func,is_linear_model,optim_opt,arg_zero_solver)
if nargin<11
    arg_zero_solver=1;
end
debug=false;
switch is_linear_model
    case 1
% %         ys0=0*ys0;
        % compute the constant
        [resid,PD]=resid_func(ys0,x_ss,[],param,def,[],[]);
        if any(resid)
            ys=ys0-pinv(full(PD))*resid; % <-- ys=ys0-PD\const;
            % The generalized inverse works better when the system is
            % singular as it is the case when the model is solved with
            % unit roots. But it might be slower and so I don't know
            % whether I should put a switch...
        else
            ys=ys0;
        end
        residuals=resid_func(ys,x_ss,[],param,def,[],[]);
        exitflag=max(abs(residuals))<=optim_opt.TolFun;
        if ~exitflag
            if debug
                disp('using lsqnonlin or fsolve on a linear model') %#ok<UNRCH>
            end
            [ys,exitflag]=steady_state_solver_engine(ys0,debug);
        end
    otherwise
        % Here is why we need good initial values. You cannot, say start at
        % zero for a variable in logs...
        residuals=resid_func(ys0,x_ss,[],param,def,[],[]);
        if isnan(residuals)
            exitflag=-4;
            ys=nan;
        elseif max(abs(residuals))<=optim_opt.TolFun
            ys=ys0;
            exitflag=1;
        else
            [ys,exitflag]=steady_state_solver_engine(ys0,debug);
        end
end
if ~any(isnan(ys)) && exitflag==1
    % maybe I should not do the following? especially here?
    ys(abs(ys)<1e-12)=0;
    retcode=0;
else
    retcode=1; % steady state cannot be solved
end

    function [x1,exitflag]=steady_state_solver_engine(x0,debug)
        if debug
            optim_opt.Display='iter';
        else
            optim_opt.Display='none';
        end
        switch arg_zero_solver
            case 1
                % call to lsqnonlin
                %------------------
                [x1,resnorm,residuals,exitflag]=lsqnonlin(resid_func,x0,[],[],optim_opt,...
                    x_ss,[],param,def,[],[]);  %#ok<*ASGLU>
            case 2
                % call to fsolve
                %------------------
                [x1,fval,exitflag]=fsolve(resid_func,x0,optim_opt,...
                    x_ss,[],param,def,[],[]);
                resnorm=norm(fval);
            otherwise
                error('arg_zero_solver must be either 1 or 2')
        end
        if ismember(exitflag,[1:4,-3]) && ~any(isnan(x1))
            if resnorm > sqrt(optim_opt.TolFun)
                exitflag=inf;
            else
                exitflag=1;
            end
        end
    end
end

