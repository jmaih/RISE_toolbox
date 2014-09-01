function [x1,f1,H,issue,viol,final_obj,final_funevals]=...
    estimation_wrapper(obj,action,x0,lb,ub,funevals)
persistent general_restrictions linear_restricts ngen_restr nobj ...
    estim_blocks lb_short ub_short npar_short oldobj old_funevals
if isempty(action)
    action='eval';
end
nonlcon=obj(1).routines.nonlinear_restrictions;
nconst=obj(1).number_of_restrictions.nonlinear;
if isempty(linear_restricts)
    if isempty(obj(1).linear_restrictions_data)
        % the model has not been estimated for the posterior mode yet. In
        % that case load the restrictions
        %-----------------------------------------------------------------
        obj=setup_linear_restrictions(obj);
        obj=setup_general_restrictions(obj);
    end
    % In order to avoid recomputing the stationarity status we
    % create a copy of the original object and use it throughout
    oldobj=obj;
    old_funevals=funevals;
    linear_restricts=oldobj(1).linear_restrictions_data;
    nobj=numel(oldobj);
    general_restrictions=cell(1,nobj);
    for iii=1:nobj
        general_restrictions{iii}=oldobj(iii).general_restrictions_data;
    end
    estim_blocks=oldobj(1).options.estim_blocks;
    if ~isempty(estim_blocks)
        estim_blocks=create_estimation_blocks(oldobj(1),estim_blocks);
    end
    lb_short=linear_restricts.a2tilde_func(lb);
    ub_short=linear_restricts.a2tilde_func(ub);
    npar_short=size(lb_short,1);
    bad=lb_short>ub_short;
    if any(bad)
        tmp=lb_short;
        lb_short(bad)=ub_short(bad);
        ub_short(bad)=tmp(bad);
    end
end
violLast=[];
xLast=[];
viol=[];
if isempty(x0)
    x0=lb_short+(ub_short-lb_short).*rand(npar_short,1);
else
    x0=linear_restricts.a2tilde_func(x0);
end
switch action
    case 'estimate'
        PROBLEM_=struct('objective',@fh_wrapper,...
            'x0',x0,...
            'lb',lb_short,...
            'ub',ub_short,...
            'nonlcon',@nonlcon_with_gradient,... % the nonlinear constraints restrictions take the same inputs as fh_wrapper
            'options',oldobj(1).options.optimset,...
            'solver',oldobj(1).options.optimizer);
        
        [x1,f1,H,issue]=optimization.estimation_engine(PROBLEM_,...
            oldobj(1).options.hessian_type,estim_blocks);
    case 'eval'
        [f1,retcode_0]=fh_wrapper(x0);
        x1=x0;
        H=[];
        issue=retcode_0;
        viol=nonlcon_with_gradient(x0);
    otherwise
        error(['unknown type of action ',action])
end
% extend the output
%-------------------
x1 = linear_restricts.a_func(x1);
if ~isempty(H)
    % the second arguments indicates that it is a covariance term
    H = linear_restricts.a_func(H,true);
end
final_obj=obj;
final_funevals=old_funevals;

    function [minus_log_post,retcode]=fh_wrapper(xtilde)
        % this function returns the minimax if there are many objects
        % evaluated simultaneously
        
        % all objects are evaluated at the same point
        % the reason you want to output the object here is because it potentially
        % contains crucial information going forward. In particular, it contains
        % information about whether the model is stationary or not.
        
        % expand x before doing anything
        %-------------------------------
        x=linear_restricts.a_func(xtilde);
        
        fval=oldobj(1).options.estim_penalty*ones(1,nobj);
        for mo=1:nobj
            [fval(mo),~,~,~,retcode,oldobj(mo)]=log_posterior_kernel(oldobj(mo),x);
            if retcode
                break
            end
        end
        old_funevals=old_funevals+1;
        % Now take the negative for minimization
        minus_log_post=-min(fval);
        % nonlinear constraints might be incompatible with blockwise
        % optimization. In that case, it is better to compute the
        % restriction violation while evaluating the objective and save
        % the results to pass on to the optimizer when it calls the
        % nonlinear constraints.
        if retcode
            if isempty(ngen_restr)
                % then we still have not found a good starting point
                % and therefore, there is no point in looking at the
                % constraints.
            else
                violLast=ones(nconst+ngen_restr,1)*realmax/(nconst+ngen_restr);
            end
        else
            violLast=mynonlinear_constraints(x,oldobj,true);
        end
        % update this element right here, so that it is ready when
        % calling nonlcon_with_gradient
        xLast=x;
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
            viol_general=[];
            viol_simple=[];
            for iobj=1:nobj
                vv=nonlcon(obj(iobj).parameter_values);
                viol_simple=[viol_simple,vv(:)]; %#ok<AGROW>
                if ~isempty(general_restrictions{iobj})
                    vv=general_restrictions{iobj}(obj(iobj));
                    viol_general=[viol_general,vv(:)]; %#ok<AGROW>
                end
            end
            if isempty(ngen_restr)
                ngen_restr=numel(viol_general);
            end
            viol=[viol_simple(:);viol_general(:)];
            xLast=x;
            violLast=viol;
        end
    end

    function [viol,grad]=nonlcon_with_gradient(xtilde)
        grad=[];
        % expand x before doing anything
        %-------------------------------
        x=linear_restricts.a_func(xtilde);
        viol=mynonlinear_constraints(x,oldobj);
    end
end
