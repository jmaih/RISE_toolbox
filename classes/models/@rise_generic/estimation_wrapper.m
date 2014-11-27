function [x1,f1,H,issue,viol,obj,funevals]=...
    estimation_wrapper(obj,action,x0,lb,ub,funevals)
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
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

if isempty(action)
    action='eval';
end
nonlcon=[];
nconst=0;
if isfield(obj(1).routines,'nonlinear_restrictions')
    nonlcon=obj(1).routines.nonlinear_restrictions;
    nconst=obj(1).number_of_restrictions.nonlinear;
end

if isempty(obj(1).linear_restrictions_data)
    % the model has not been estimated for the posterior mode yet. In
    % that case load the restrictions
    %-----------------------------------------------------------------
    obj=setup_linear_restrictions(obj);
    obj=setup_general_restrictions(obj);
end

linear_restricts=obj(1).linear_restrictions_data;
nobj=numel(obj);
ngenrest=0;
c=cell(1,nobj);
for iii=1:nobj
    % check the filtering level required under estimation
    %-----------------------------------------------------
    if ~isempty(obj(iii).general_restrictions_data)
        ngenrest=ngenrest+1;
        req=obj(iii).general_restrictions_data();
        obj.options.kf_filtering_level=req.kf_filtering_level;
        c{iii}=obj(iii).options.estim_penalty_factor;
    end
end

estim_blocks=obj(1).options.estim_blocks;
if ~isempty(estim_blocks)
    estim_blocks=create_estimation_blocks(obj(1),estim_blocks);
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
            'options',obj(1).options.optimset,...
            'solver',obj(1).options.optimizer);
        
        [x1,f1,H,issue]=optimization.estimation_engine(PROBLEM_,...
            obj(1).options.hessian_type,estim_blocks);
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
        
        fval=obj(1).options.estim_penalty*ones(1,nobj);
        
        % general restrictions are sometimes infeasible. Rather than
        % imposing a barrier that will be difficult to cross if the initial
        % constraints are not satified, we use a penalty function approach.
        for mo=1:nobj
            [fval(mo),~,~,~,retcode,obj(mo)]=log_posterior_kernel(obj(mo),x);
            if retcode
                break
            end
        end
        
        if ~retcode && ngenrest
            % add a penalty for the restrictions violation
            %----------------------------------------------
            g=evaluate_nonlinear_restrictions(obj);
            for mo=1:nobj
                if ~isempty(g{mo})
                    fval(mo)=fval(mo)+utils.estim.penalize_violations(g{mo},c{mo});
                end
            end
        end
        funevals=funevals+1;
        % Now take the negative for minimization
        minus_log_post=-min(fval);
        % nonlinear constraints might be incompatible with blockwise
        % optimization. In that case, it is better to compute the
        % restriction violation while evaluating the objective and save
        % the results to pass on to the optimizer when it calls the
        % nonlinear constraints.
        if retcode
            violLast=ones(nconst,1)*realmax/(nconst);
            % update this element right here, so that it is ready when
            % calling nonlcon_with_gradient.
            xLast=x;
        else
            violLast=mynonlinear_constraints(x,obj,true);
        end
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
        % expand x before doing anything
        %-------------------------------
        x=linear_restricts.a_func(xtilde);
        viol=mynonlinear_constraints(x,obj);
    end
end
