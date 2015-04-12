function [ys,retcode]=solve_steady_state(ys0,def,pp,resid_func,...
    optim_opt,arg_zero_solver)

% solve_steady_state solves the steady state
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
% More About
% ------------
%
% Examples
% ---------
%
% See also:

if nargin<6
    arg_zero_solver=1;
end
debug=optim_opt.debug;
TolFun=optim_opt.TolFun;
if TolFun<1e-6
    TolFun=1e-6;
end

nregs=size(ys0,2);
ys=ys0;
retcode=0;
%{
% compute the unique steady state based on the ergodic distribution
%------------------------------------------------------------------
if optim_opt.is_unique
    % if the steady state is unique and the probabilities are endogenous,
    % then the steady state transition matrix cannot be computed
    % independently of the value of the endogenous variables. This is
    % something to correct as soon as possible
    [TransMat,retcode]=compute_steady_state_transition_matrix(...
        optim_opt.trans_mat_func,...
        ys0(:,1),pp(:,1),def{1},optim_opt.exo_nbr);
    if ~retcode
        [pp_unique,def_unique,retcode]=...
            dsge_tools.ergodic_parameters(TransMat.Qinit,def,pp);
        % override the parameters and the definitions as they might be used
        % for further processing in case the steady state is imposed
        %------------------------------------------------------------------
        pp=pp_unique(:,ones(1,nregs));
        def=cellfun(@(x)def_unique,def,'uniformOutput',false);
    end
end
%}
if ~retcode
    exitflag=1;
    for ireg=1:nregs
        % the parameters may have been overriden in case the steady state
        % is unique
        pp_i=pp(:,ireg);
        def_i=def{ireg};
        if exitflag==1
            if optim_opt.is_linear_model
                % compute the constant
                [resid,PD]=resid_func(ys0(:,ireg),pp_i,def_i);
                if any(resid)
                    ys(:,ireg)=ys0(:,ireg)-pinv(full(PD))*resid; % <-- ys=ys0-PD\const;
                    % The generalized inverse works better when the system is
                    % singular as it is the case when the model is solved with
                    % unit roots. But it might be slower and so I don't know
                    % whether I should put a switch...
                end
                residuals=resid_func(ys(:,ireg),pp_i,def_i);
                exitflag=max(abs(residuals))<=TolFun;
            else
                % Here is why we need good initial values. You cannot, say start at
                % zero for a variable in logs...
                residuals=resid_func(ys0(:,ireg),pp_i,def_i);
                nan_residuals=isnan(residuals);
                imag_residuals=logical(imag(residuals));
                if any(nan_residuals|imag_residuals)
                    if debug
                        if any(nan_residuals)
                            culprits=find(nan_residuals);
                            type='nan';
                        else
                            type='imaginary';
                            culprits=find(imag_residuals);
                        end
                        disp(culprits)
                        disp(['the equations # above have ',type,' residuals in the computation of the steady state'])
                        keyboard
                    end
                    exitflag=-4;
                elseif max(abs(residuals))<=TolFun
                    exitflag=1;
                else
                    [ys(:,ireg),exitflag]=steady_state_solver_engine(ys0(:,ireg));
                end
            end
        end
        if optim_opt.is_unique && nregs>1
            ys=ys(:,ones(1,nregs));
            break
        end
    end
    
    if isreal(y) && ~any(isnan(ys(:))) && all(isfinite(ys(:))) && exitflag==1
        % maybe I should not do the following? especially here?
        ys(abs(ys)<1e-12)=0;
        retcode=0;
    else
        retcode=1; % steady state cannot be solved
    end
end

    function [x1,exitflag]=steady_state_solver_engine(x0)
        if debug
            optim_opt.Display='iter';
        else
            optim_opt.Display='none';
        end
        switch arg_zero_solver
            case 1
                % call to lsqnonlin
                %------------------
                [x1,resnorm,residuals,exitflag]=lsqnonlin(resid_func,x0(:),[],[],optim_opt,pp_i,def_i);  %#ok<*ASGLU>
            case 2
                % call to fsolve
                %------------------
                [x1,fval,exitflag]=fsolve(resid_func,x0(:),optim_opt,pp_i,def_i);
                resnorm=norm(fval);
            otherwise
                error('arg_zero_solver must be either 1 or 2')
        end
        exitflag=utils.optim.exitflag(exitflag,x1,resnorm,sqrt(TolFun));
        if exitflag==1
            x1=reshape(x1,size(x0));
        end
    end
end

