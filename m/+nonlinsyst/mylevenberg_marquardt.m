function [x,fval,exitflag] = mylevenberg_marquardt(f,x,opts,varargin)

warningstate1 = []; warningstate2 = []; warningstate3 = [];

sqrtEps=sqrt(eps);

lambda=0.01;

fcount=0;

nv=length(x);

I=eye(nv);

objfun=@objective;

opts=nonlinsyst.set_defaults(opts);

% differencing interval for numerical gradient
%---------------------------------------------
delta = opts.TolX;

tvec=delta*I;

is_verbose=strcmp(opts.Display,'iter');

[fval,sumSq,Jacobian,badJacobian] = objfun(x);

is_success=~badJacobian;

if ~is_success
    
    error('bad initial jacobian...')
    
end

% last successful Jacobian: but what if the first jacobian is bad?
%-----------------------------------------------------------------
oldJacob=Jacobian;

[dx,~,gradF]=set_step(x,fval,Jacobian,is_success);

iter=0;

% Initial 1st order optimality with safeguard to prevent value of 0, used in stopping tests
relFactor = max(norm(gradF,Inf),sqrtEps);

[done,exitflag,iter] = check_convergence(gradF,relFactor,opts,iter,badJacobian, ...
    x,[],sqrtEps,dx,[],fcount);

% first iteration information
if is_verbose
    formatstr          = ' %5.0f       %5.0f   %13.6g    %12.3g %12.6g   %12.6g\n';
    fprintf( ...
        ['\n                                        First-Order                    Norm of \n', ...
        ' Iteration  Func-count    Residual       optimality      Lambda           step\n']);
    fprintf(formatstr,iter,fcount,sumSq,norm(gradF,Inf),lambda,nan);
end

while ~done
    
    [dx,xnew]=set_step(x,fval,Jacobian,is_success);
    
    % If the previous step wasn't successful, we don't need to evaluate the
    % Jacobian until we're sure that the latest step is a good one. Only
    % evaluate the cost function in that case.
    if ~is_success
        
        [fvalnew,trialSumSq] = objfun(xnew);
        
    else
        
        [fvalnew,trialSumSq,Jacobian,badJacobian] = objfun(xnew);
        
    end
        
    goodTrial = isfinite(trialSumSq) && isreal(trialSumSq);
    
    if goodTrial && trialSumSq < sumSq
        % We've reduced the sum squared error and the function vector is
        % defined at the new trial point.
        fval = fvalnew;         
        
        x = xnew;
        
        if is_success
            
            % Reduce LM parameter, only if previous step was good
            lambda = 0.1*lambda;    
            
        end
        % If Jacobian and previous step(s) unsuccessful, recompute the
        % Jacobian right here, right now
        if ~is_success
            % avoid recomputing the reference guy...
            [~,~,Jacobian,badJacobian] = objfun(xnew,fvalnew);
            
        end
        
        is_success = true;        % Successful step, new Jacobian needs to be computed
        
        gradF = Jacobian'*fval;
        
        % Print iterative display
        if is_verbose
            
            fprintf(formatstr,iter,fcount,trialSumSq,norm(gradF,Inf),lambda,norm(dx));
            
        end
        
        % Check Termination Criteria
        [done,exitflag,iter] = check_convergence(gradF,relFactor,opts,iter,...
            badJacobian,x,trialSumSq,sqrtEps,dx,sumSq,fcount);
        
        sumSq = trialSumSq;         % Update sum of squares after convergence test
        
    else
        
        lambda = 10*lambda;              % Increase LM parameter
        
        is_success = false;          % Unsuccessful step, no need to re-compute the Jacobian
        
        % The step is too small, cannot proceed or too many function evals
        if norm(dx) < opts.TolX*(sqrtEps + norm(xnew))
            
            exitflag = 4;
            
            done = true;
            
        elseif fcount > opts.MaxFunEvals
            
            exitflag = 0;
            
            done = true;
            
        end
        
    end
    
end

    function [dx,xnew,gradF]=set_step(x,fval,Jacobian,is_success)
        
        % Compute LM step
        if is_success
            
            oldJacob=Jacobian;
            
        end
        
        gradF = Jacobian'*fval;
        
        control_warnings(true)
                
        %===============================
        JJ=oldJacob.'*oldJacob;
                
        dx=-(JJ + lambda*I)\(oldJacob.'*fval);
        %===============================
        
        control_warnings(false)
        
        xnew = x + dx;                  % Update X with computed step
        
    end

    function [o,sumSquares,g,badJacobian]=objective(x,o)
        
        if nargin<2
            
            o=f(x,varargin{:});
            
            fcount=fcount+1;
            
        end
        
        sumSquares=nonlinsyst.norm(o);
        
        if nargout>2
            
            g = o(:,ones(nv,1));
            
            f0=o;
            
            for i=1:nv
                
                f1=objective(x+tvec(:,i));
                
                g(:,i) = (f1-f0)/delta;
                
            end
            
            badJacobian=~all(isfinite(g(:)));
            
        end
        
    end

    function control_warnings(do_set)
        
        if do_set
            % Disable the warnings about conditioning for singular and
            % nearly singular matrices
            
            warningstate1 = warning('off','MATLAB:nearlySingularMatrix');
            warningstate2 = warning('off','MATLAB:singularMatrix');
            warningstate3 = warning('off','MATLAB:rankDeficientMatrix');
            
        else
            
            % Restore the warning states to their original settings
            warning(warningstate1)
            warning(warningstate2)
            warning(warningstate3)
            
        end
        
    end

end

function [done,exitflag,iter] = check_convergence(gradF,relFactor,opts,...
    iter,badJacobian,x,newF,sqrtEps,step,oldF,fcount)

% Initialize these quantities in case no criteria are met.
done = false;

exitflag = [];

% tolOpt: tolerance used when checking the 1st-order optimality
tolOpt = 1e-4 * opts.TolFun;

if badJacobian

    exitflag = 2;
    
    done = true;

elseif norm(gradF,Inf) < tolOpt * relFactor

    exitflag = 1;
    
    done = true;

elseif iter > 0

    if norm(step) < opts.TolX*(sqrtEps + norm(x))
    
        exitflag = 4;
        
        done = true;
    
    elseif abs(newF - oldF) <= opts.TolFun*oldF
    
        exitflag = 3;
        
        done = true;
    
    elseif (fcount > opts.MaxFunEvals)||(iter >=opts.MaxIter)
    
        exitflag = 0;
        
        done = true;
    
    else
        
        iter = iter + 1;
    
    end
    
else
    
    iter = iter + 1;

end

end
