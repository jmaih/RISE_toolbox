function [obj,structural_matrices,retcode]=compute_steady_state(obj,varargin)
% compute_steady_state - computes the steady state of a dsge model
%
% Syntax
% -------
% ::
%
%   [obj,structural_matrices,retcode]=compute_steady_state(obj,varargin)
%
% Inputs
% -------
%
% - **obj** [rise|dsge]: model file
%
% - **varargin** []: usual optional arguments
%
% Outputs
% --------
%
% - **obj** [rise|dsge]: model file
%
% - **structural_matrices** [struct]: structure containing various
%   important elements for the solution of the system
%
% - **retcode** [numeric]: 0 if there was problem computing the steady
%   state.
%
% More About
% ------------
%
% - There are 2 cases to consider:
%   - The user does not provide any steady state equations: RISE will
%   attempt to solve the steady state using a vector of zeros as initial
%   guess. It makes life easy if the user provides the status of the
%   variables in the system i.e. whether they grow linear or log-linearly.
%   - The user provide some equations for solving the steady state. This is
%   done in two ways:
%       1-) the steady_state_model block: the variables that do not appear
%           in the block will be initialized at zero. Some parameters can
%           also be computed inside the block. The user can define an
%           optimization to solve for a subset of steady state values
%           inside the block. The block has three attributes
%           (i) imposed(default=false): RISE computes the solution at the
%           specified point without checking that the point solves for the
%           steady state
%           (ii) unique (default=false): RISE computes the steady state at
%           the ergodic distribution of the parameters. In case the
%           probabilities are endogenous, the ergodic distribution of the
%           parameters is itself a function of the steady state of the
%           variables.
%           (iii) loop(default=false): RISE considers the equations
%           calculating the steady state as true and just solves for the
%           missing variables by looping over the steady state program. The
%           user can then use the values pushed into the steady state
%           program to calculate the steady state for the included
%           variables.
%       2-) the steady state file: The user writes a function which can be
%           called in two possible ways: (i) [vnames,info]=ssfile(); In
%           this case the first output argument is the list of variables
%           for which the user computes the steady state; the second output
%           is a structure with fields unique, imposed and loop
%           just as in the case of the steady state model. (ii) The other
%           call to the function is [y,newp,retcode]=ssfile(y,p,d,id,obj).
%           In this case, the first input (y) is the vector of steady
%           states, which is updated and returned as the first output. The
%           locations of the modifications are indicated by the fourth
%           input (id), which is computed based on the list of the
%           variables in vnames above. As for the other outputs, p is a
%           structure with parameters, d is a structure with definitions,
%           obj is the model object in case the user needs some further
%           information for computing the steady state. In case some
%           parameters are computed in the steady state file, they should
%           be returned in the structure "newp". The last output "retcode"
%           indicates whether no problem was encoutered in the computation
%           of the steady state (retcode=0) or the opposite (retcode =
%           any number different from 0).
%
% - Writing the steady state file in this ways makes it possible to use the
%   same whether there is regime switching or not.
%
% Examples
% ---------
%
% See also:

if isempty(obj)
    
    obj=struct('steady_state_file','',...
        'steady_state_use_steady_state_model',true,...
        'steady_state_solver','lsqnonlin',...
        'steady_state_algorithm',{{'levenberg-marquardt',2*.005}});

    return
    
end

structural_matrices=struct();

if ~isempty(varargin)
    
    obj=set(obj,varargin{:});
    
end

[obj,retcode]=compute_definitions(obj);

if retcode
    
    return
    
end

[obj,sscode,blocks]=prepare_steady_state_program(obj);

number_of_regimes=obj.markov_chains.small_markov_chain_info.regimes_number;

endo_nbr=obj.endogenous.number;

exo_nbr=sum(obj.exogenous.number);

unsolved=true(1,endo_nbr);

x=zeros(exo_nbr,1);

if isempty(obj.old_solution)
    
    g0=zeros(endo_nbr,1);
    
    g0(obj.endogenous.is_log_var,:)=1;
    
    y0=g0;
    
    y0=y0(:,ones(1,number_of_regimes));
    
    g0=g0(:,ones(1,number_of_regimes));
    
else
    
    g0=cell2mat(obj.old_solution.bgp);
    
    y0=cell2mat(obj.old_solution.ss);
    
end

p=obj.parameter_values;

d=obj.solution.definitions;

r=nan(endo_nbr,number_of_regimes);

y=y0;

g=g0;

if obj.is_unique_steady_state
    
    converged=false;
    
    while ~converged
        
        [p_unic,def_unic,TransMat,retcode]=ergodic_parameters(y0(:,1));
        
        if retcode
            
            break
            
        end
        
        [y(:,1),g(:,1),~,p_unic,retcode]=run_one_regime(y0(:,1),g0(:,1),p_unic,x,...
            def_unic,sscode,unsolved,blocks);
       
        converged=retcode||max(abs(y0(:,1)-y(:,1)))<=obj.options.fix_point_TolFun;
        
        y0(:,1)=y(:,1);
        
    end
    
    if ~retcode
        
        y=y(:,ones(1,number_of_regimes));
        
        g=g(:,ones(1,number_of_regimes));
        
        % update the parameters since we are not going into the steady
        % state for every regime
        %--------------------------------------------------------------
        update_changed_parameters(p_unic);
        
        % recompute residuals for all regimes
        %-------------------------------------
        for i2=1:number_of_regimes
            
            r(:,i2)=blocks.residcode(y(:,i2),g(:,i2),x,p(:,i2),d{i2});
            
        end
        
    end
    
else
    
    for istate=1:number_of_regimes
        
        [y(:,istate),g(:,istate),r(:,istate),p(:,istate),retcode]=...
            run_one_regime(y0(:,istate),g0(:,istate),p(:,istate),x,...
            d{istate},sscode,unsolved,blocks);

        if retcode
            
            break
            
        end
        
        if ~isempty(sscode) && sscode.is_imposed_steady_state
            
            % update the parameters since we are not going into the steady
            % state for every regime
            %--------------------------------------------------------------
            update_changed_parameters(p(:,1));
            
            y=y(:,ones(1,number_of_regimes));
            
            g=g(:,ones(1,number_of_regimes));
            
            % recompute residuals for the remaining regimes
            %-----------------------------------------------
            for i2=2:number_of_regimes
                
                r(:,i2)=blocks.residcode(y(:,i2),g(:,i2),x,p(:,i2),d{i2});
                
            end
            
            break
            
        end
        
    end
    
    if ~retcode
        
        [TransMat,retcode]=compute_steady_state_transition_matrix(...
            obj.routines.transition_matrix,y(:,1),p(:,1),d{1},...
            sum(obj.exogenous.number));
        
        
    end
    
end


if ~retcode
    
    obj.parameter_values=p;
    
    obj.solution.transition_matrices=TransMat;
    
    
    % prepare output
    %-----------------
    structural_matrices.user_resids=r;
    
    structural_matrices.transition_matrices=obj.solution.transition_matrices;
    
    y=sparse(y);
    
    g=sparse(g);
    
    ss_tvp=y(obj.endogenous.is_affect_trans_probs,:);
    
    bad=any(abs(bsxfun(@minus,ss_tvp,ss_tvp(:,1)))>1e-9,2);
    
    if any(bad)
        
        bad_endo_vars=get(obj,'endo_list(affect_trans_probs)');
        
        bad_endo_vars=bad_endo_vars(bad);
        
        disp(bad_endo_vars)
        
        error(['The variables above affect the transition probabilities but ',...
            'do not have the same steady state in each regime'])
        
    end
    
    obj.solution.ss=mat2cell(y,endo_nbr,ones(1,number_of_regimes));
    
    obj.solution.bgp=mat2cell(g,endo_nbr,ones(1,number_of_regimes));
    
end

    function update_changed_parameters(ptarg)
        
        if ~isempty(sscode) && any(sscode.is_param_changed) && sscode.is_unique_steady_state
            % Push the parameters around if they have changed in the steady
            % state programs
            %--------------------------------------------------------------
            changelocs=find(sscode.is_param_changed);
            
            p(changelocs,2:number_of_regimes)=ptarg(changelocs)*ones(1,number_of_regimes-1);
            
        end
        
    end

    function [pp_unique,def_unique,TransMat,retcode]=ergodic_parameters(y)
        
        [TransMat,retcode]=compute_steady_state_transition_matrix(...
            obj.routines.transition_matrix,y,p(:,1),...
            d{1},sum(obj.exogenous.number));
        
        if retcode
            
            pp_unique=[];
            
            def_unique=[];
            
        else
            
            [pp_unique,def_unique,retcode]=...
                dsge_tools.ergodic_parameters(TransMat.Qinit,d,p);
            
        end
        
    end

end



function flag=check_residuals(r,tol)

flag=all(isfinite(r)) && max(abs(r))<tol;

end



function [y,g,r,p,retcode]=run_one_regime(y0,g0,p,x,d,sscode,unsolved,dq_blocks)

is_log_var=dq_blocks.is_log_var;

solve_bgp=dq_blocks.solve_bgp;

auxcode=dq_blocks.auxcode;

optimopt=dq_blocks.optimopt;

if ~isempty(sscode)
    if sscode.is_loop_steady_state
        
        [y,g,r,retcode]=loop_over(y0,g0,sscode.subset);
        
    else
        
        [y,g,p,r,retcode]=grand_sscode(y0,g0,p,d);
        
        if sscode.is_imposed_steady_state || check_residuals(r,optimopt.TolFun)
            
            return
            
        else
            
            if retcode
                
                unsolved(:)=true;
                
                [y,g,r,retcode]=divide_and_conquer(y,g,p,x,d,unsolved,...
                    dq_blocks);
                
            else
                
                if any(unsolved)
                    
                    [y,g,r,retcode]=divide_and_conquer(y,g,p,x,d,unsolved,...
                        dq_blocks);
                    
                end
                
            end
            
        end
        
    end
    
else
    
    [y,g,r,retcode]=divide_and_conquer(y0,g0,p,x,d,unsolved,dq_blocks);
    
end



    function [y,g,r,retcode]=loop_over(y0,g0,subset)
        
        if islogical(subset)
            
            nv=sum(subset);
            
        else
            
            nv=numel(subset);
            
        end
        
        arg_zero_solver=dq_blocks.arg_zero_solver;
        
        % pay some respect to log vars
        %-----------------------------
        log_status=is_log_var(subset);
        
        ssc0=y0(subset);
        
        if solve_bgp
            
            ssc0=[ssc0;g0(subset)];
            
            log_status=[log_status,log_status]; 
            
        end
        
        ssc0(log_status) = log(ssc0(log_status));
        
        [ssc1,retcode]=optimization_center(@concentrated_residuals,ssc0,...
            arg_zero_solver,optimopt);
        
        [r,y,g]=concentrated_residuals(ssc1);
        
        function [r,y,g]=concentrated_residuals(ssc)
            
            ssc(log_status) = exp(ssc(log_status));
            
            ytmp=y0; 
            
            ytmp(subset)=ssc(1:nv);
            
            gtmp=g0;
            
            if solve_bgp
                
                gtmp(subset)=ssc(nv+1:end);
                
            end
            
            [y,g,p,r,retcode]=grand_sscode(ytmp,gtmp,p,d);
            
            if ~all(isfinite(r))
                
                retcode=1;
                
                r(~isfinite(r))=1;
                
            end
            
        end
        
    end



    function [y,g,p,r,retcode]=grand_sscode(y0,g0,p,d)
        
        yg0=y0;
        
        g = g0;
        
        if solve_bgp
            
            yg0=yg0+g0*1i;
            
        end
        
        [yg,p,retcode]=sscode.func(yg0,p,d);
        
        if ~retcode
            
            [yg,retcode]=auxcode(yg,x,p,d);
            
        end
        
        y=real(yg);
        
        if solve_bgp
            
            g=imag(yg);
            
        end
        
        if retcode
            
            r=ones(size(y));
            
        else
            
            r=dq_blocks.residcode(y,g,x,p,d);
            
        end
        
    end

end



function [y,g,r,retcode]=divide_and_conquer(y,g,p,x,d,unsolved,blocks)

optimopt=blocks.optimopt;

% try a quick exit
%-----------------
r=blocks.residcode(y,g,x,p,d);

if check_residuals(r,optimopt.TolFun)
    
    retcode=0;
    
    return
    
end

nblocks=numel(blocks.variables);

arg_zero_solver=blocks.arg_zero_solver;

is_log_var=blocks.is_log_var;

solve_bgp=blocks.solve_bgp;

fixed_vars=find(~unsolved);

sseqtns=blocks.sseqtns;

debug=blocks.debug;

shift=blocks.solve_bgp_shift;

fast=nblocks==numel(sseqtns);

for iblock=1:nblocks
    
    varblk=blocks.variables{iblock};
    
    fixed=ismember(varblk,fixed_vars);
    
    varblk(fixed)=[];
    
    nv=numel(varblk);
    
    if nv==0
        continue
    end
    
    eqtnblk=blocks.equations{iblock};
    
    neqtn=numel(eqtnblk);
    
    log_status=is_log_var(varblk);
    
    if solve_bgp
        
        log_status=[log_status,log_status]; %#ok<AGROW>
        
    end
    
    z0=y(varblk);
    
    if solve_bgp
        
        z0=[z0;g(varblk)]; %#ok<AGROW>
        
    end
    
    z0(log_status)=log(z0(log_status));
    
    [z,retcode]=optimization_center(@small_system,z0,arg_zero_solver,optimopt);

    z(log_status)=exp(z(log_status));
    
    z(abs(z) <= optimopt.TolX) = 0;
    
    y(varblk)=z(1:nv);
    
    if solve_bgp
        
        g(varblk)=z(nv+1:end);
        
    end
    
    if debug
        
        blocks.endo_list(varblk)
        
        blocks.dynamic(eqtnblk)
        
        disp([blocks.endo_list(varblk)',num2cell(full([y(varblk),g(varblk)]))])
        
    end
    
    if retcode
        
        return
        
    end
    
end

r=blocks.residcode(y,g,x,p,d);

    function rs=small_system(zyg)
        
        zyg(log_status)=exp(zyg(log_status));
        
        y(varblk)=zyg(1:nv);
        
        if solve_bgp
            
            g(varblk)=zyg(nv+1:end);
            
        end
        
        if fast
            
            rs=sseqtns{iblock}(y,g,x,p,d);
            
        else
            
            rs=zeros(neqtn,1);
            
            for ieqtn=1:neqtn
                
                rs(ieqtn)=sseqtns{eqtnblk(ieqtn)}(y,g,x,p,d);
                
            end
            
        end
        
        if solve_bgp
            
            yshift=y;
            
            yshift(~is_log_var)=y(~is_log_var)+shift*g(~is_log_var);
            
            yshift(is_log_var)=y(is_log_var).*g(is_log_var).^shift;
            
            if fast
                
                rs2=sseqtns{iblock}(yshift,g,x,p,d);
                
            else
                
                rs2=rs;
                
                for ieqtn=1:neqtn
                    
                    rs2(ieqtn)=sseqtns{eqtnblk(ieqtn)}(yshift,g,x,p,d);
                    
                end
                
            end
            
            rs=[rs;rs2];
            
        end
        
%         if any(~isfinite(rs)),rs=10000*ones(size(rs));end
        
    end

end



function [ssc1,retcode]=optimization_center(objfun,ssc0,arg_zero_solver,optimopt)

if strcmp(arg_zero_solver,'lsqnonlin')
    % call to lsqnonlin
    %------------------
    [ssc1,resnorm,~,exitflag]=lsqnonlin(objfun,ssc0,[],[],optimopt);
    
elseif strcmp(arg_zero_solver,'fsolve')
    % call to fsolve
    %------------------
    [ssc1,fval,exitflag]=fsolve(objfun,ssc0,optimopt);
    
    resnorm=norm(fval);
    
elseif strcmp(arg_zero_solver,'fminunc')
    % call to fsolve
    %------------------
    [ssc1,resnorm,exitflag]=fminunc(@minimizer,ssc0,optimopt);
    
else
    
    [arg_zero_solver,vargs]=utils.code.user_function_to_rise_function(arg_zero_solver);
    
    [ssc1,resnorm,exitflag]=arg_zero_solver(objfun,ssc0,optimopt,vargs{:});
    
end
% N.B: I use sqrt(tol) below on purpose!!!
%-----------------------------------------
sqrt_tol=optimopt.TolFun;

exitflag=utils.optim.exitflag(exitflag,ssc1,resnorm,sqrt_tol);
%--------------------------------------------------------------
retcode=1-(exitflag==1);

    function fval=minimizer(ssc)
        
        fval=objfun(ssc);
        
        fval=norm(fval);
        
    end

end



function [obj,sscode,blocks]=prepare_steady_state_program(obj)
% steady state file
%------------------
sscode=struct();

sscode.is_param_changed=false(1,sum(obj.parameters.number));

if ~isempty(obj.options.steady_state_file) && ...
        ~isa(obj.options.steady_state_file,'function_handle')
    
    if ~ischar(obj.options.steady_state_file)
        
        error('steady state file must be a function handle or a char')
        
    end
    
    obj.options.steady_state_file=str2func(obj.options.steady_state_file);
    
end

sscode=[];

default_sstate_attributes=struct('unique',false,'imposed',false,...
    'loop',false);

if isempty(obj.steady_state_2_model_communication)
    
    if ~isempty(obj.options.steady_state_file)
        
        sscode=struct();
        
        % memo will follow
        %------------------
        sscode.func=obj.options.steady_state_file;
        
        % call the function with the object only
        %----------------------------------------
        [var_names,updated_param_list,new_settings]=sscode.func(obj);
        
        isUpdate=false(1,obj.parameters.number(1));
        if ~isempty(updated_param_list)
            locs=locate_variables(updated_param_list,obj.parameters.name);
            isUpdate(locs)=true;
        end
        
        sscode.is_param_changed=isUpdate;
        
        ff=fieldnames(new_settings);
        
        for ifield=1:numel(ff)
            
            if ~isfield(default_sstate_attributes,ff{ifield})
                error(['The expected fields of a steady state file are : ',...
                    '"unique","imposed" and "loop"'])
            end
            
            default_sstate_attributes.(ff{ifield})=new_settings.(ff{ifield});
            
        end
        
        sscode.is_imposed_steady_state=default_sstate_attributes.imposed;
        
        sscode.is_unique_steady_state=default_sstate_attributes.unique;
        
        sscode.is_loop_steady_state=default_sstate_attributes.loop;
        
        user_endo_ids=locate_variables(var_names,obj.endogenous.name);
        
        sscode.solved=ismember(obj.endogenous.name,var_names);
        
        sscode.ids=user_endo_ids;
        
        % the following two lines could consume a lot of memory!!!
        % If the problem becomes accute, one will just have to drop
        % steady_state_2_model_communication and rebuild the structure
        % for each new parameter vector.
        %--------------------------------------------------------
        sscode.func=memoize_steady_statefile(obj,...
            sscode.func,sscode.ids,obj.parameters.name,...
            obj.definitions.name);
        
    end
    
    % steady state model
    %--------------------
    if isempty(sscode)
        
        if ~isempty(obj.routines.steady_state_model) && ...
                ~isempty(obj.routines.steady_state_model.code) && ...
                obj.options.steady_state_use_steady_state_model
            
            sscode=struct();
            
            sscode.func=memoize_steady_statemodel(obj.routines.steady_state_model);
            
            sscode.is_param_changed=obj.is_param_changed_in_ssmodel;
            
            sscode.solved=obj.auxiliary_variables.ssmodel_solved;
            
            sscode.is_imposed_steady_state=obj.is_imposed_steady_state;
            
            sscode.is_unique_steady_state=obj.is_unique_steady_state;
            
            sscode.is_loop_steady_state=obj.is_loop_steady_state;
            
        end
        
    end
    
    if ~isempty(sscode)
        
        % after having parsed the ssfile and the ssmodel, we can do this
        %----------------------------------------------------------------
        
        errmsg=['With an initial guess for the steady state, ',...
            'the parameters cannot be modified in the steady state ',...
            'model or file'];
        
        if ~(sscode.is_imposed_steady_state||...
                sscode.is_loop_steady_state) && ...
                any(sscode.is_param_changed)
            
            error(errmsg)
            
        end
        
        if  sscode.is_loop_steady_state
            % subset of variables to loop over
            %---------------------------------
            sscode.subset=~sscode.solved & ~obj.endogenous.is_auxiliary;
            
        end
        
        
        obj.steady_state_2_model_communication=sscode;
        
    end
    
else
    
    sscode=obj.steady_state_2_model_communication;
    
end


if isempty(obj.steady_state_2_blocks_optimization)
    
    % optimization setup
    %--------------------
    optimopt=obj.options.optimset;
    
    optimopt.TolX=1e-12;
    
    optimopt.TolFun=1e-12;
    
    optimopt.Algorithm=obj.options.steady_state_algorithm;
    
    if obj.options.debug
        
        optimopt.Display='iter';
        
    else
        
        optimopt.Display='none';
        
    end
    
    % auxiliary variables
    %--------------------
    planner_routines=[];
    
    if obj.is_optimal_policy_model
        
        planner_routines=struct(...
            'planner_static_mult_support',obj.routines.planner_static_mult_support,...
            'planner_static_mult',obj.routines.planner_static_mult);
        
    end
    
    auxcode=auxiliary_memoizer(...
        obj.routines.shadow_steady_state_auxiliary_eqtns,...
        planner_routines,...
        obj.endogenous.is_lagrange_multiplier);
    
    blocks=struct('variables',{obj.steady_state_blocks.variables},...
        'is_log_var',obj.endogenous.is_log_var,...
        'solve_bgp',obj.options.solve_bgp,...
        'solve_bgp_shift',obj.options.solve_bgp_shift,...
        'sseqtns',{obj.routines.static.separate_blocks},...
        'debug',obj.options.debug,...
        'equations',{obj.steady_state_blocks.equations},...
        'endo_list',{obj.endogenous.name},...
        'dynamic',{obj.equations.dynamic},...
        'residcode',obj.routines.static.one_block,...
        'auxcode',auxcode,...
        'arg_zero_solver',obj.options.steady_state_solver,...
        'optimopt',optimopt);
    
    obj.steady_state_2_blocks_optimization=blocks;
    
else
    
    blocks=obj.steady_state_2_blocks_optimization;
    
end


end



function memo=memoize_steady_statemodel(ssmodelcode)

memo=@memo_engine;

    function varargout=memo_engine(y,p,d)
        % the steady state model returns the steady state and the
        % parameters but not the retcode, unlike the steady state file. So
        % it has to be constructed.
        nout=nargout;
        
        varargout=cell(1,nout);
        
        ss=y;
        
        x=[]; s0 =[]; s1 =[];
        
        [varargout{1:2}]=utils.code.evaluate_functions(ssmodelcode,y,x,ss,p,d,s0,s1);
        
        retcode=0;
        
        if any(isnan(varargout{1}))||any(isinf(varargout{1}))
            retcode=1;
        end
        
        varargout{3}=retcode;
    end

end



function memo=memoize_steady_statefile(obj,ssfilecode,id,pnames,defnames)

if ischar(ssfilecode)
    ssfilecode=str2func(ssfilecode);
end

memo=@engine;

    function varargout=engine(y,p,d)
        
        [varargout{1:nargout}]=ssfilecode(obj,y,...
            vector2struct(pnames,p),...
            vector2struct(defnames,d),id);
     
        % the second output includes changed parameters and it has to be
        % dealt with consequently. In particular, we want to make it so
        % that there is no difference between ssfile and ssmodel
        %-----------------------------------------------------------------
        newp=varargout{2};
        
        if ~isempty(newp)
            
            fields=fieldnames(newp);
            
            locs=locate_variables(fields,pnames);
            
            for ifield=1:numel(locs)
                
                p(locs(ifield))=newp.(fields{ifield});
                
            end
            
        end
        
        varargout{2}=p;
    end

    function pnew=vector2struct(vnames,v)
        pnew=struct();
        for ip=1:numel(v)
            pnew.(vnames{ip})=v(ip);
        end
    end

end



function auxcode=auxiliary_memoizer(aux_ssfunc,planner_routines,is_lagr_mult)

auxcode=@auxiliary_evaluation;

    function [yg,retcode]=auxiliary_evaluation(yg,x,p,d)
        
        y=real(yg);
        
        g=imag(yg);
        
        retcode=0;
        
        if ~isempty(aux_ssfunc) && (isa(aux_ssfunc,'function_handle')||...
                (isstruct(aux_ssfunc)&&~isempty(aux_ssfunc.code)))

            y=utils.code.evaluate_functions(aux_ssfunc,y,g,x,p,d);
            
            g=utils.code.evaluate_functions(aux_ssfunc,g,g,x,p,d);
            
        end
        
        if ~isempty(planner_routines)
            
            % TODO: Think of a way to compute the growth rate of the
            % lagrange multipliers besides doing it numerically through
            % divide and conquer.
            
            siz=planner_routines.planner_static_mult_support{1};
            
            pos=planner_routines.planner_static_mult_support{2};
            
            vals=zeros(size(pos));
            
            vals(pos)=utils.code.evaluate_functions(planner_routines.planner_static_mult,...
                y,x,ss,p,d,[],[]); % y, x, ss, param, def, s0, s1
            
            vals=reshape(vals,siz);
            
            wx=vals(end,:);
            
            vals(end,:)=[];
            
            if any(wx)
                
                multvals=vals.'\wx.';
                
            else
                
                nmults=sum(is_lagr_mult);
                
                multvals=zeros(nmults,1);
                
            end
            
            y(is_lagr_mult)=multvals;
            
        end
        
        yg=y+g*1i;
        
    end

end