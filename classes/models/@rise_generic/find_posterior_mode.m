function [x1,f1,H,x0,f0,viol,funevals,issue,obj]=find_posterior_mode(obj,x0,lb,ub,...
    basics,general_restrictions,gen_restr_args)
nobj=numel(obj);
if nobj==0
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    x1=struct();
    return
end

optim_options=obj(1).options.optimset;
optimizer=obj(1).options.optimizer;
hessian_type=obj(1).options.hessian_type;
estim_blocks=obj(1).options.estim_blocks;
if ~isempty(estim_blocks)
    estim_blocks=create_estimation_blocks(obj(1),estim_blocks);
end

[nonlcon,nconst]=reprocess_nonlinear_restrictions(obj(1).parameter_restrictions);
% the functions in obj(1).parameter_random_draws_restrictions only check
% whether the restrictions are violated or not but do not give the strength
% of the violation, which we need in order to apply DEB. So we need to use
% the reprocessed form.

% initialize the number of function calls
%----------------------------------------
funevals=0;

Nsim=max(1,obj(1).options.estim_parallel);
% now shorten everything
%-----------------------
x0=basics.a2tilde_func(x0);
npar_short=size(x0,1);
x0=[x0,nan(npar_short,Nsim-1)];
f0=nan(1,Nsim);
lb=basics.a2tilde_func(lb);
ub=basics.a2tilde_func(ub);
bad=lb>ub;
if any(bad)
    tmp=lb;
    lb(bad)=ub(bad);
    ub(bad)=tmp(bad);
end
% all objects are evaluated at the same point. Without a second argument,
% this is exactly what will happen.
[~,f0(1),~,retcode0,viol]=big_wrapper(x0(:,1));
if retcode0||any(viol>0)
    % first check constraint violations
    f0(1)=obj(1).options.estim_penalty;
end
% the objects are automatically updated and potentially contain crucial
% information going forward. In particular, they contain information
% about whether the models are stationary or not.

if f0(1)<obj(1).options.estim_penalty
    beg=2;
else
    beg=1;
end

%% disable those elements
warning('off','MATLAB:nearlySingularMatrix')
warning('off','MATLAB:illConditionedMatrix')

fprintf(1,'%s\n','Looking for good enough start values. Please wait...');
for ii=beg:Nsim
    NotDone=true;
    iter=0;
    while NotDone
        [xtest,ftest,retcode]=utils.estim.generate_starting_point(@big_wrapper);
        if ftest<obj(1).options.estim_penalty
            NotDone=false;
            f0(ii)=ftest;
            x0(:,ii)=xtest;
        end
        iter=iter+1;
        if iter>=obj(1).options.estim_max_trials
            error([mfilename,':: No admissible starting value after ',...
                int2str(obj(1).options.estim_max_trials),' trials'])
        else
            fprintf(1,'%3.0d :: %s\n',iter,utils.error.decipher(retcode));
        end
    end
    disp(['Starting value # ',int2str(ii),' found after ',int2str(iter),' iterations'])
    ratio=ii/Nsim;
    fprintf(1,'%s\n',['...', num2str(100*ratio),'% done']);
end

[x1,f1,H,issue,viol,obj]=big_wrapper(x0,'estimate');
% both vector x1 and H here are short: extend them
%-------------------------------------------------
x1 = basics.a_func(x1);
H = basics.a_func(H,true); % the second arguments indicates that it is a covariance term
% also resize x0, which can be a matrix
%--------------------------------------
x0 = basics.a_func(x0);

viol=viol(viol>0);

warning('on','MATLAB:nearlySingularMatrix')
warning('on','MATLAB:illConditionedMatrix')

    function [x1,f1,H,issue,viol,finalobj]=big_wrapper(x0,action)
        if nargin<2
            action='eval';
            if nargin<1
                x0=lb+(ub-lb).*rand(npar_short,1);
            end
        end
        violLast=[];
        xLast=[];
        viol=[];
        ngen_restr=[];
        switch action
            case 'estimate'
                PROBLEM_=struct('objective',@fh_wrapper,...
                    'x0',x0,...
                    'lb',lb,...
                    'ub',ub,...
                    'nonlcon',@nonlcon_with_gradient,... % the nonlinear constraints restrictions take the same inputs as fh_wrapper
                    'options',optim_options,...
                    'solver',optimizer);
                
                [x1,f1,H,issue]=optimization.estimation_engine(PROBLEM_,hessian_type,estim_blocks);
            case 'eval'
                [f1,retcode_0]=fh_wrapper(x0);
                x1=x0;
                H=[];
                issue=retcode_0;
                viol=nonlcon_with_gradient(x0);
            otherwise
                error(['unknown type of action ',action])
        end
        finalobj=obj;
        
        function [minus_log_post,retcode]=fh_wrapper(xtilde)
            % this function returns the minimax if there are many objects
            % evaluated simultaneously
            
            % all objects are evaluated at the same point
            % the reason you want to output the object here is because it potentially
            % contains crucial information going forward. In particular, it contains
            % information about whether the model is stationary or not.
            
            % expand x before doing anything
            %-------------------------------
            x=basics.a_func(xtilde);
            
            fval=obj(1).options.estim_penalty*ones(1,nobj);
            for mo=1:nobj
                [fval(mo),~,~,~,retcode,obj(mo)]=log_posterior_kernel(obj(mo),x);
                if retcode
                    break
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
                if isempty(ngen_restr)
                    % then we still have not found a good starting point
                    % and therefore, there is no point in looking at the
                    % constraints.
                else
                    violLast=ones(nconst+ngen_restr,1)*realmax/(nconst+ngen_restr);
                end
            else
                violLast=mynonlinear_constraints(x,obj);
            end
            % update this element right here, so that it is ready when
            % calling nonlcon_with_gradient
            xLast=x;
        end
        
        function viol=mynonlinear_constraints(x,obj) %#ok<INUSL>
            % this function assumes that the first argument has already
            % been pushed into the second one... I will need to make sure
            % that this is actually the case.
            viol_general=[];
            viol_simple=[];
            for iobj=1:nobj
                vv=nonlcon(obj(iobj).parameter_values);
                viol_simple=[viol_simple,vv(:)]; %#ok<AGROW>
                if ~isempty(general_restrictions{iobj})
                    vv=general_restrictions{iobj}(obj(iobj),...
                        gen_restr_args{iobj}{:});
                    viol_general=[viol_general,vv(:)]; %#ok<AGROW>
                end
            end
            if isempty(ngen_restr)
                ngen_restr=numel(viol_general);
            end
            viol=[viol_simple(:);viol_general(:)];
        end
        
        function [viol,grad]=nonlcon_with_gradient(xtilde)
            grad=[];
            % expand x before doing anything
            %-------------------------------
            x=basics.a_func(xtilde);
            if isequal(x,xLast)
                viol=violLast;
            else
                % % % % % % % % % % % %                 warning('evaluating constraints before computing log-posterior ')
                thisobj=assign_estimates(obj,x);
                % % % % % % % % % % % %                 thisobj=filter(thisobj);
                viol=mynonlinear_constraints(x,thisobj);
            end
        end
    end

end

function [nonlcon,nconst]=reprocess_nonlinear_restrictions(nonlcon)
% transform the nonlinear constraints. I would like to keep the
% flexibility of knowning what parameters enter the constraints and
% so I do not do this in format parameters
nconst=numel(nonlcon);
for iconstr=1:nconst
    % remove semicolon
    nonlcon{iconstr}=strrep(nonlcon{iconstr},';','');
    % now remove inequalities
    cutoff_type={'>=','<=','>','<','='};
    for itype=1:numel(cutoff_type)
        cutoff_locs=strfind(nonlcon{iconstr},cutoff_type{itype});
        if ~isempty(cutoff_locs)
            cutoff_type=cutoff_type{itype};
            break
        end
    end
    if ~isempty(cutoff_locs)
        span=length(cutoff_type);
        left=nonlcon{iconstr}(1:cutoff_locs-1);
        right=nonlcon{iconstr}(cutoff_locs+span:end);
        switch cutoff_type
            case '>='
                nonlcon{iconstr}=[right,'-(',left,')-eps;'];
            case '<='
                nonlcon{iconstr}=[left,'-(',right,')-eps;'];
            case '>'
                nonlcon{iconstr}=[right,'-(',left,');'];
            case '<'
                nonlcon{iconstr}=[left,'-(',right,');'];
            case '='
                nonlcon{iconstr}=['abs(',left,'-(',right,'))-eps;'];
        end
    end
end
nonlcon=cell2mat(nonlcon(:)');
nonlcon=nonlcon(1:end-1);
if isempty(nonlcon)
    nonlcon='0';
    nconst=1;
end
nonlcon=str2func(['@(M)[',nonlcon,']']);
end
% function rst=push_time_series(C)
% rst=ts(C{1},C{2},C{3});
% end
