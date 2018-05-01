function [x1,f1,H,issue,viol,obj,funevals]=...
estimation_wrapper(obj,action,x0,lb,ub,funevals)

% estimation_wrapper -- environment for the estimation of parameters
%
% ::
%
%
%   [x1,f1,H,issue,viol,obj,funevals]=estimation_wrapper(obj,action,x0,lb,ub,funevals)
%
% Args:
%
%    - **obj** [rise|dsge|svar|rfvar]: model object
%
%    - **action** ['draw'|{'eval'}|'estimate']: intended action
%
%    - **x0** [[]|d x 1 vector]: paramter vector for start of the optimization
%    or for log-posterior-kernel evaluation
%
%    - **lb** [d x 1 vector]: lower bound of the parameters to estimate
%
%    - **ub** [d x 1 vector]: upper bound of the parameters to estimate
%
%    - **funevals** [[]|numeric]: function evaluations
%
% Returns:
%    :
%
%    - **x1** [[]|d x 1 vector]: final paramter vector
%
%    - **f1** [numeric]: value of objective function evaluated at x1
%
%    - **H** [d x d x 2 array]: Hessian, such that H(:,:,1) is the hessian
%    computed returned from the optimizer and H(:,:,2) is the hessian computed
%    numerically by RISE.
%
%    - **issue** [char|cellstr]: list of problems encountered during the
%    process
%
%    - **viol** [[]|vector]: restriction violations
%
%    - **obj** [rise|dsge|svar|rfvar]: model object possibly modified
%
%    - **funevals** [[]|numeric]: function evaluations
%
% Note:
%
% Example:
%
%    See also:

if isempty(action)
    
    action='eval';
    
end

if strcmp(action,'draw')
    
    npar=numel(lb);
    
    x0=lb+(ub-lb).*rand(npar,1);
    
    action='eval';
    
end

nonlcon=[];

nconst=0;

if isfield(obj(1).routines,'nonlinear_restrictions')
    
    nonlcon=obj(1).routines.nonlinear_restrictions;
    
    nconst=obj(1).number_of_restrictions.nonlinear;
    
end

% the model has not been estimated for the posterior mode yet. In
% that case load the restrictions
%-----------------------------------------------------------------
obj=setup_restrictions(obj);

nobj=numel(obj);

ngenrest=0;

c=cell(1,nobj);

for iii=1:nobj
    
    if ~isempty(obj(iii).general_restrictions_data)
        
        ngenrest=ngenrest+1;
        
        c{iii}=obj(iii).options.estim_penalty_factor;
        
    end
    
end

estim_penalty = -abs(obj(1).options.estim_penalty);

estim_barrier=obj(1).options.estim_barrier;

estim_blocks=obj(1).options.estim_blocks;

if ~isempty(estim_blocks)
    
    estim_blocks=create_estimation_blocks(obj(1),estim_blocks);
    
end

violLast=[];

xLast=[];

viol=[];

switch action
    
    case 'estimate'
        
        PROBLEM_=struct('objective',@fh_wrapper,...
            'x0',x0,...
            'lb',lb,...
            'ub',ub,...
            'nonlcon',@nonlcon_with_gradient,... % the nonlinear constraints restrictions take the same inputs as fh_wrapper
            'options',obj(1).options.optimset,...
            'solver',obj(1).options.optimizer);
        
        [x1,f1,H,issue]=optimization.estimation_engine(PROBLEM_,estim_blocks);
        
    case 'eval'
        
        [f1,retcode_0,max_pen]=fh_wrapper(x0);
        
        x1=x0;
        
        H=[];
        
        issue=retcode_0;
        
        viol=nonlcon_with_gradient(x0);
        
        if isempty(viol)
            
            viol=max_pen;
            
        else
            
            viol=viol+max_pen;
            
        end
        
    otherwise
        
        error(['unknown type of action ',action])
        
end

    function [minus_log_post,retcode,max_pen]=fh_wrapper(xtilde)
        % this function returns the minimax if there are many objects
        % evaluated simultaneously
        
        % all objects are evaluated at the same point
        % the reason you want to output the object here is because it potentially
        % contains crucial information going forward. In particular, it contains
        % information about whether the model is stationary or not.
        
        % untransform in the presence of dirichlet right here and right
        % now. Otherwise there will be a wrong prior evaluation and further
        % problems down the road
        %------------------------------------------------------------------
        x=unstransform_parameters(obj(1),xtilde);
        
        fval=estim_penalty*ones(1,nobj);
        
        is_processed=false;
        
        if estim_barrier
            
            violLast=mynonlinear_constraints(x,obj,false);
            
            if any(violLast>0)
                
                is_processed=true;
                
            end
            
        end
        
        max_pen=0;
        
        retcode=0;
        
        if ~is_processed
            % general restrictions are sometimes infeasible. Rather than
            % imposing a barrier that will be difficult to cross if the initial
            % constraints are not satified, we use a penalty function approach.
            for mo=1:nobj
                
                [fval(mo),loglik,~,~,retcode,obj(mo)]=log_posterior_kernel(obj(mo),x);
                
                if retcode
                    
                    break
                    
                end
                
                if obj.is_fisher
                    % return the likelihood only
                    fval(mo)=loglik;
                    
                end
                
            end
            
        end
        
        % nonlinear constraints might be incompatible with blockwise
        % optimization. In that case, it is better to compute the
        % restriction violation while evaluating the objective and save
        % the results to pass on to the optimizer when it calls the
        % nonlinear constraints.
        if retcode
            
            violLast=ones(nconst,1)*realmax/(nconst);
            
        else
            
            if ngenrest
                % add a penalty for the restrictions violation
                %----------------------------------------------
                g=evaluate_general_restrictions(obj);
                
                for mo=1:nobj
                    
                    if ~isempty(g{mo})
                        
                        if estim_barrier && any(g{mo}>0)
                            
                            fval(:)=estim_penalty;
                            
                            is_processed=true;
                            
                        end
                        
                        if ~is_processed
                            
                            this_pen=utils.estim.penalize_violations(g{mo},max(abs(fval(mo)),c{mo}));
                            % violations of the restrictions decrease the value
                            % of the objective function
                            fval(mo)=fval(mo)-this_pen;
                            
                            max_pen=max(max_pen,this_pen);
                            
                        end
                        
                    end
                    
                end
                
            end
            
            if ~is_processed
                
                violLast=mynonlinear_constraints(x,obj,true);
                
            end
            
        end
        
        % make xLast ready for nonlcon_with_gradient.
        %----------------------------------------------
        xLast=x;
        
        % Now take the negative for minimization
        %----------------------------------------
        minus_log_post=-min(fval);
        
        % update the number of function evaluations
        %-------------------------------------------
        funevals=funevals+1;
        
    end

    function viol=mynonlinear_constraints(x,obj,params_pushed)
        
        if nargin<3
            
            params_pushed=false;
            
        end
        % this function assumes that the first argument has already
        % been pushed into the second one... I will need to make sure
        % that this is actually the case.
        if isequal(x,xLast)
            
            viol=violLast;
            
        else
            
            if ~params_pushed
                
                obj=assign_estimates(obj,x);
                
            end
            
            viol_simple=[];
            
            for iobj=1:nobj
                
                if ~isempty(nonlcon)
                    
                    vv=nonlcon(obj(iobj).parameter_values);
                    
                    viol_simple=[viol_simple,vv(:)]; %#ok<AGROW>
                    
                end
                
            end
            
            viol=viol_simple(:);
            
            xLast=x;
            
            violLast=viol;
            
        end
        
    end

    function [viol,grad]=nonlcon_with_gradient(xtilde)
        
        grad=[];
        
        % unstransform xtilde into x before doing anything
        %--------------------------------------------------
        x=unstransform_parameters(obj(1),xtilde);
    
        viol=mynonlinear_constraints(x,obj);

    end

end
