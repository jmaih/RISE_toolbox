function [obj,structural_matrices,retcode]=compute_steady_state(obj,varargin)
% compute_steady_state - computes the steady state of a dsge model
%
% ::
%
%
%   [obj,structural_matrices,retcode]=compute_steady_state(obj,varargin)
%
% Args:
%
%    - **obj** [rise|dsge]: model file
%
%    - **varargin** []: usual optional arguments
%
% Returns:
%    :
%
%    - **obj** [rise|dsge]: model file
%
%    - **structural_matrices** [struct]: structure containing various
%      important elements for the solution of the system
%
%    - **retcode** [numeric]: 0 if there was problem computing the steady
%      state.
%
% Note:
%
%    - There are 2 cases to consider:
%      - The user does not provide any steady state equations: RISE will
%      attempt to solve the steady state using a vector of zeros as initial
%      guess. It makes life easy if the user provides the status of the
%      variables in the system i.e. whether they grow linear or log-linearly.
%      - The user provide some equations for solving the steady state. This is
%      done in two ways:
%          1-) the steady_state_model block: the variables that do not appear
%              in the block will be initialized at zero. Some parameters can
%              also be computed inside the block. The user can define an
%              optimization to solve for a subset of steady state values
%              inside the block. The block has three attributes
%              (i) imposed(default=false): RISE computes the solution at the
%              specified point without checking that the point solves for the
%              steady state
%              (ii) unique (default=false): RISE computes the steady state at
%              the ergodic distribution of the parameters. In case the
%              probabilities are endogenous, the ergodic distribution of the
%              parameters is itself a function of the steady state of the
%              variables.
%              (iii) loop(default=false): RISE considers the equations
%              calculating the steady state as true and just solves for the
%              missing variables by looping over the steady state program. The
%              user can then use the values pushed into the steady state
%              program to calculate the steady state for the included
%              variables.
%          2-) the steady state file: The user writes a function which can be
%              called in two possible ways: (i) [vnames,info]=ssfile(); In
%              this case the first output argument is the list of variables
%              for which the user computes the steady state; the second output
%              is a structure with fields unique, imposed and loop
%              just as in the case of the steady state model. (ii) The other
%              call to the function is [y,newp,retcode]=ssfile(y,p,d,id,obj).
%              In this case, the first input (y) is the vector of steady
%              states, which is updated and returned as the first output. The
%              locations of the modifications are indicated by the fourth
%              input (id), which is computed based on the list of the
%              variables in vnames above. As for the other outputs, p is a
%              structure with parameters, d is a structure with definitions,
%              obj is the model object in case the user needs some further
%              information for computing the steady state. In case some
%              parameters are computed in the steady state file, they should
%              be returned in the structure "newp". The last output "retcode"
%              indicates whether no problem was encoutered in the computation
%              of the steady state (retcode=0) or the opposite (retcode =
%              any number different from 0).
%
%    - Writing the steady state file in this ways makes it possible to use the
%      same whether there is regime switching or not.
%
% Example:
%
%    See also:

% possible fields are "imposed", "unique", "loop".
% The default value for all of those is false.
%   - "imposed": This tells RISE not to check that this is actually solves
%       the steady state. Hence, RISE will attempt to approximate the model
%       around the chosen point
%   - "unique": This tells RISE that the steady state is the same across
%       all regimes. RISE will call the function only once but instead of
%       using just any parameter vector, it will use the ergodic mean of
%       the parameter vector (mean across all regimes).
%   - "loop": This tells RISE that if the user was given the steady state
%       values for some of the variables, he would be able to compute the
%       steady state for the remaining variables. RISE will then exploit
%       this information to optimize over the variables that the user needs
%       for computing the steady state.

% - As it seems, the true jacobian only works in linear models. In
% nonlinear models, somehow, the finite difference approach yields better
% results.

if isempty(obj)
    
    obj=cell(0,4);
    
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

[obj,sscode,blocks,unsolved]=prepare_steady_state_program(obj);

number_of_regimes=obj.markov_chains.small_markov_chain_info.regimes_number;

endo_nbr=obj.endogenous.number;

exo_nbr=sum(obj.exogenous.number);

x=zeros(exo_nbr,1);

if isempty(obj.old_solution)
    
    g0=zeros(endo_nbr,1);
    
    g0(obj.endogenous.is_log_var,:)=1;
    
    y0=g0;
    
    y0=y0(:,ones(1,number_of_regimes));
    
    g0=g0(:,ones(1,number_of_regimes));
    
else
    
    big2small=true;
    
    map_=map_regimes(obj,big2small);
    
    g0=cell2mat(obj.old_solution.bgp);
    
    g0=g0(:,map_);
    
    y0=cell2mat(obj.old_solution.ss);
    
    y0=y0(:,map_);
    
end

p=obj.parameter_values;

d=obj.solution.definitions;

r=nan(endo_nbr,number_of_regimes);

y=y0;

g=g0;

blocks.y0=y0;

blocks.g0=g0;

swapping_func=blocks.swapping_func;

extend_param_func=blocks.extend_param_func;

if obj.options.steady_state_unique
    
    converged=false;
    
    iloop=0;
        
    p_start=p(:,1);
    
    def_start=d{1};
    
    TransMat=[];
    
    while ~converged
        
        iloop=iloop+1;
        
        [y(:,1),g(:,1),~,p_start,retcode]=run_one_regime(y0(:,1),g0(:,1),p_start,x,...
            def_start,sscode,unsolved,blocks);
        
        % exit if there is a problem (retcode>0) or if we have converged
        %----------------------------------------------------------------
        converged=retcode||max(abs(y0(:,1)-y(:,1)))<=obj.options.fix_point_TolFun;
        
        if ~retcode
            
            [p_unic,def_unic,TransMat,retcode]=ergodic_parameters(y(:,1));
            
            if retcode
                
                break
                
            end
            
            y0(:,1)=y(:,1);
            
            g0(:,1)=g(:,1);
            
            p_start=p_unic;
            
            def_start=def_unic;
            
        end
        
    end
    
    if ~retcode
        
        % swap for the output
        %---------------------
        [y(:,1),p_unic]=swapping_func(y(:,1),p_unic);
        
        % expand
        %--------
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
        
        % swap for the output
        %---------------------
        [y(:,istate),p(:,istate)]=swapping_func(y(:,istate),p(:,istate));
        
        if istate==1
            
            [TransMat,retcode]=compute_steady_state_transition_matrix(...
                obj.routines.transition_matrix,y(:,1),p(:,1),d{1},...
                sum(obj.exogenous.number));
            
        end
        
    end
    
end


if ~retcode
    
    obj.parameter_values=p;
    
    obj.solution.transition_matrices=TransMat;
    
    
    % prepare output
    %-----------------
    structural_matrices.user_resids=r;
    
    if ~isreal(r)||any(isnan(r(:)))
        
        retcode = 11;
        
        return
        
    end
    
    structural_matrices.transition_matrices=obj.solution.transition_matrices;
    
    y=sparse(y);
    
    g=sparse(g);
    
    ss_tvp=y(obj.endogenous.is_affect_trans_probs,:);
    
    bad=any(abs(bsxfun(@minus,ss_tvp,ss_tvp(:,1)))>1e-9,2);
    
    if any(bad) && obj.options.debug
        
        bad_endo_vars=get(obj,'endo_list(affect_trans_probs)');
        
        bad_endo_vars=bad_endo_vars(bad);
        
        disp(bad_endo_vars)
        
        warning(['The variables above affect the transition probabilities but ',...
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
            
            p(changelocs,:)=ptarg(changelocs)*ones(1,number_of_regimes);
            
        end
        
        % alternatively, apply the swapping...
        %---------------------------------------
        p=extend_param_func(p,ptarg);
        
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
    
    swapping_func=dq_blocks.swapping_func;
    
    if sscode.is_loop_steady_state && ~dq_blocks.solve_linear
        % no point in doing loops if the model is linear
        
        [y,g,r,retcode]=loop_over(y0,g0,sscode.subset);
                
    else
        
        [y,g,p,r,retcode]=grand_sscode(y0,g0,p,d);
        
        if sscode.is_imposed_steady_state || check_residuals(r,optimopt.TolFun)
            % do nothing
        else
            
            if retcode
                
                unsolved(:)=true;
                
                [y,g,p,r,retcode]=divide_and_conquer(y,g,p,x,d,unsolved,...
                    dq_blocks);
                
            else
                
                if any(unsolved)||~all(isfinite(y))||~all(isfinite(g))
                    
                    [y,g,p,r,retcode]=divide_and_conquer(y,g,p,x,d,unsolved,...
                        dq_blocks);
                    
                end
                
            end
            
        end
        
    end
    
else
    
    [y,g,p,r,retcode]=divide_and_conquer(y0,g0,p,x,d,unsolved,dq_blocks);
    
end

if solve_bgp
    
    r=r(1:numel(y));
    
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
        
        % try a quick exit
        %-----------------
        r=evaluate_all_residuals(dq_blocks.residcode,y0,g0,x,p,d,swapping_func);
        
        if check_residuals(r,optimopt.TolFun)
            
            can_proceed=true;
            
            if dq_blocks.solve_bgp
                
                yshift=create_shift(y0,g0,is_log_var,dq_blocks.solve_bgp_shift);
                
                rs2=evaluate_all_residuals(dq_blocks.residcode,yshift,g0,x,p,d,swapping_func);
                
                r=[r;rs2];
                
                can_proceed=check_residuals(rs2,optimopt.TolFun);
                
            end
            
            if can_proceed
                
                y=y0;
                
                g=g0;
                
                retcode=0;
                
                % possibly update the parameter
                %------------------------------
                if any(sscode.is_param_changed)
                    
                    [~,p,retcode]=sscode.func(y,p,d);
                    
                end
                
                if retcode==0
                    
                    return
                    
                end
                
            end
            
        end
        
        yg0=y0;
        
        g = g0;
        
        if solve_bgp
            
            yg0=[yg0,g0];
            
        end
        
        [yg,p,retcode]=sscode.func(yg0,p,d);
        
        if ~retcode
            
            [yg,retcode]=auxcode(yg,x,p,d);
            
        end
        
        y=yg(:,1);
        
        if solve_bgp
            
            g=yg(:,2);
            
        end
        
        if retcode
            
            nyg=(1+dq_blocks.solve_bgp)*size(y,1);
            
            r=ones(nyg,1);
            
        else
            
            r=evaluate_all_residuals(dq_blocks.residcode,y,g,x,p,d,swapping_func);
            
            if dq_blocks.solve_bgp
                
                yshift=create_shift(y,g,is_log_var,dq_blocks.solve_bgp_shift);
                
                rs2=evaluate_all_residuals(dq_blocks.residcode,yshift,g,x,p,d,swapping_func);
                
                r=[r;rs2];
                
            end
            
        end
        
    end

end



function [y,g,p,r,retcode]=divide_and_conquer(y,g,p,x,d,unsolved,blocks)

optimopt=blocks.optimopt;

swapping_func=blocks.swapping_func;

solve_bgp=blocks.solve_bgp;

% try a quick exit
%-----------------
r=evaluate_all_residuals(blocks.residcode,y,g,x,p,d,swapping_func);

shift=blocks.solve_bgp_shift;

is_log_var=blocks.is_log_var;

if check_residuals(r,optimopt.TolFun)
    
    can_proceed=true;
    
    if solve_bgp
        
        yshift=create_shift(y,g,is_log_var,shift);
        
        rs2=evaluate_all_residuals(blocks.residcode,yshift,g,x,p,d,swapping_func);
        
        can_proceed=check_residuals(rs2,optimopt.TolFun);
        
        r=[r;rs2];
        
    end
    
    if can_proceed
        
        retcode=0;
        
        return
        
    end
    
end

sseqtns=blocks.sseqtns;

ssjacobians=blocks.ssjacobians;

debug=blocks.debug;

if debug
    
    y0=y;
    
    g0=g;
    
end

if blocks.solve_linear
    % exploit jacobian information and exit quickly
    %-----------------------------------------------
    if any(is_log_var)
        
        error('log variables not permitted in linear models')
        
    end
    
    iblock=1;
    
    [r,retcode]=solve_linearly();
    
    return
    
end

nblocks=numel(blocks.variables);

arg_zero_solver=blocks.arg_zero_solver;

fixed_vars=find(~unsolved);

for iblock=1:nblocks
    
    [orig_varlk,varblk,nv,fixed_level,fixed_growth]=load_relevant();
    
    if nv>0
        
        level_vars=varblk(~fixed_level);
        
        nvl=numel(level_vars);
        
        growth_vars=varblk(~fixed_growth);
        
        log_status=is_log_var(level_vars);
        
        if solve_bgp
            
            log_status=[log_status,...
                is_log_var(growth_vars)]; %#ok<AGROW>
            
        end
        
        z0=y(level_vars);
        
        if solve_bgp
            
            z0=[z0;g(growth_vars)]; %#ok<AGROW>
            
        end
        
        z0(log_status)=log(z0(log_status));
        
        [z,retcode]=optimization_center(@small_system,z0,arg_zero_solver,optimopt);
        
        z(log_status)=exp(z(log_status));
        
        try
            z=round(z,12);
        catch
            z(abs(z) <= optimopt.TolX) = 0;
        end
        
        y(level_vars)=z(1:nvl);
        
        if solve_bgp
            
            g(growth_vars)=z(nvl+1:end);
            
        end
        
    else
        
        retcode=0;
        
    end
    
    if debug
        
        eqtnblk=blocks.equations{iblock};
        
        blocks.endo_list(orig_varlk)
        
        blocks.dynamic(eqtnblk)
        
        disp([blocks.endo_list(orig_varlk)',...
            num2cell(full([y(orig_varlk),g(orig_varlk),y0(orig_varlk),g0(orig_varlk)]))])
        
    end
    
    if retcode
        
        return
        
    end
    
end

r=evaluate_all_residuals(blocks.residcode,y,g,x,p,d,swapping_func);

    function [orig_varblk,varblk,nv,fixed_level,fixed_growth]=load_relevant()
        
        varblk=blocks.variables{iblock};
        
        orig_varblk=varblk;
        
        fixed=ismember(varblk,fixed_vars);
        
        fixed=fixed(:);
        
        bad_y=isnan(y(varblk));
        
        bad_g=isnan(g(varblk));
        
        fixed_level=fixed & ~bad_y;
        
        fixed_growth=fixed & ~bad_g;
        
        y(varblk(bad_y))=blocks.y0(varblk(bad_y));
        
        g(varblk(bad_g))=blocks.g0(varblk(bad_g));
        
        totally_fixed=fixed_level & fixed_growth;
        
        varblk(totally_fixed)=[];
        
        fixed_level(totally_fixed)=[];
        
        fixed_growth(totally_fixed)=[];
        
        nv=numel(varblk);
    end

    function [rupdate,retcode]=solve_linearly()
        % Since the model is linear, the Jacobian is independent of the
        % variables.
        
        J=solve_jacobian(ssjacobians{iblock},y,[],g);
        
        if solve_bgp
            
            ny=numel(y);
            
            Ay=J(:,1:ny);
            
            Ag=J(:,ny+1:end);
            
            r2=r+shift*Ay*g;
            
            yg=[y;g]-pinv(full([J;Ay,shift*Ay+Ag]))*[r;r2];
            
            yg(abs(yg)<optimopt.TolFun)=0;
            
            y=yg(1:ny);
            
            g=yg(ny+1:end);
            
        else
            
            y=J\uminus(r)+y;
            
        end
        
        rupdate=blocks.residcode(y,g,x,p,d);
        
        resnorm=norm(rupdate);
        
        exitflag=utils.optim.exitflag(1,[y;g],resnorm);
        
        retcode=1-(exitflag==1);
        
    end

    function J=solve_jacobian(func,y,yshift,g)
        
        J=utils.code.evaluate_functions(func,y,g,x,p,d);
        
        if ~isempty(yshift)
            
            J2=utils.code.evaluate_functions(func,yshift,g,x,p,d);
            
            J=[J;J2];
            
        end
        
    end

    function [rs,Jac]=small_system(zyg)
        
        is_jac=nargout>1;
        
        zyg(log_status)=exp(zyg(log_status));
        
        y(level_vars)=zyg(1:nvl);
        
        if solve_bgp
            
            g(growth_vars)=zyg(nvl+1:end);
            
        end
        
        rs=evaluate_all_residuals(sseqtns{iblock},y,g,x,p,d,swapping_func);
        
        yshift=[];
        
        if solve_bgp
            
            yshift=create_shift(y,g,is_log_var,shift);
            
            rs2=evaluate_all_residuals(sseqtns{iblock},yshift,g,x,p,d,swapping_func);
            
            rs=[rs;rs2];
            
        end
        
        if is_jac
            
            Jac=solve_jacobian(ssjacobians{iblock},y,yshift,g);
            
        end
        
        %         if any(~isfinite(rs)),rs=10000*ones(size(rs));end
        
    end

end



function yshift=create_shift(y,g,is_log_var,shift)

yshift=y;

yshift(~is_log_var)=y(~is_log_var)+shift*g(~is_log_var);

yshift(is_log_var)=y(is_log_var).*g(is_log_var).^shift;

end



function r=evaluate_all_residuals(residfunc,y,g,x,p,d,swapping_func)

if ~isempty(swapping_func)
    
    [y,p]=swapping_func(y,p);
    
end

r=residfunc(y,g,x,p,d);

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

exitflag=utils.optim.exitflag(exitflag,ssc1,resnorm);

retcode=1-(exitflag==1);

    function fval=minimizer(ssc)
        
        fval=objfun(ssc);
        
        fval=norm(fval);
        
    end

end



function [obj,sscode,blocks,unsolved]=prepare_steady_state_program(obj)

% In case of a swap, the location of the endogenous variables swapped
% remain unsolved !!!

unsolved=true(1,obj.endogenous.number);

% steady state functions (just for output)
%-----------------------------------------
if obj.warrant_setup_change
    
    obj=recreate_steady_state_functions(obj);
    
    obj.steady_state_2_blocks_optimization=[];
    
    obj.steady_state_2_model_communication=[];
    
end

% steady state file
%------------------
sscode=[];

if isempty(obj.steady_state_2_model_communication)
    
    if ~isempty(obj.options.steady_state_file)
        
        if ~isa(obj.options.steady_state_file,'function_handle')
            
            if ~ischar(obj.options.steady_state_file)
                
                error('steady state file must be a function handle or a char')
                
            end
            
            obj.options.steady_state_file=str2func(obj.options.steady_state_file);
            
        end
        
        sscode=struct();
        
        % memo will follow
        %------------------
        sscode.func=obj.options.steady_state_file;
        
        % call the function with the object only
        %----------------------------------------
        [var_names,updated_param_list]=sscode.func(obj);
        
        isUpdate=false(1,obj.parameters.number(1));
        
        if ~isempty(updated_param_list)
            
            locs=locate_variables(updated_param_list,obj.parameters.name);
            
            isUpdate(locs)=true;
            
        end
        
        sscode.is_param_changed=isUpdate;
        
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
                (isa(obj.routines.steady_state_model,'function_handle')||...
                ~isempty(obj.routines.steady_state_model.code)) && ...
                obj.options.steady_state_use_steady_state_model
            
            sscode=struct();
            
            sscode.func=memoize_steady_statemodel(obj.routines.steady_state_model);
            
            sscode.is_param_changed=obj.is_param_changed_in_ssmodel;
            
            sscode.solved=obj.auxiliary_variables.ssmodel_solved;
            
        end
        
    end
    
    if ~isempty(sscode)
        
        sscode.is_imposed_steady_state=obj.options.steady_state_imposed;
        
        sscode.is_unique_steady_state=obj.options.steady_state_unique;
        
        sscode.is_loop_steady_state=obj.options.steady_state_loop;
        
        sscode.is_fixed_steady_state=obj.options.steady_state_fixed;
        
        if sscode.is_fixed_steady_state && sscode.is_loop_steady_state
            
            error('steady state cannot be simultaneously fixed and looped')
            
        end
        
        % after having parsed the ssfile and the ssmodel, we can do this
        %----------------------------------------------------------------
        
        errmsg=['With an initial guess for the steady state, ',...
            'the parameters should not be modified in the steady state ',...
            'model or file'];
        
        if ~(sscode.is_imposed_steady_state||...
                sscode.is_loop_steady_state) && ...
                any(sscode.is_param_changed)
            
%             warning(errmsg)
            
        end
        
        if  sscode.is_loop_steady_state
            % subset of variables to loop over
            %---------------------------------
            sscode.subset=~sscode.solved & ~obj.endogenous.is_auxiliary & ...
                ~obj.endogenous.is_lagrange_multiplier;
            
        end
        
        obj.steady_state_2_model_communication=sscode;
        
    end
    
else
    
    sscode=obj.steady_state_2_model_communication;
    
end

if ~isempty(sscode) && sscode.is_fixed_steady_state
    
    unsolved(sscode.solved)=false;
    
end

if isempty(obj.steady_state_2_blocks_optimization)||...
        obj.warrant_setup_change
    
    % optimization setup
    %--------------------
    optimopt=obj.options.optimset;
    
    % we no longer overwrite the default tolerance
    %---------------------------------------------
    % optimopt.TolX=1e-12; optimopt.TolFun=1e-12;
    
    % When a solution is to be found, usually it is found very fast. When
    % there is no solution, it may take forever so we change the default
    % number of iterations from 1000 to 20
    %---------------------------------------------------------------
    %     optimopt.MaxIter=20;
    
    optimopt.Algorithm=obj.options.steady_state_algorithm;
    
    % The true jacobian is faster but it does not work as well as the
    % finite difference approximation for lsqnonlin and fsolve when the
    % model is nonlinear. In the linear case, however, we do need the
    % jacobian.
    if obj.options.steady_state_use_jacobian
        
        optimopt.Jacobian='on';
        
    else
        
        optimopt.Jacobian='off';
        
    end
    
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
            'planner_static_mult_support',{obj.routines.planner_static_mult_support},...
            'planner_static_mult',{obj.routines.planner_static_mult});
        
    end
    
    auxcode=auxiliary_memoizer(...
        obj.routines.shadow_steady_state_auxiliary_eqtns,...
        obj.routines.shadow_steady_state_auxiliary_growth_rate_eqtns,...
        planner_routines,...
        obj.endogenous.is_lagrange_multiplier);
    
    blocks=struct('variables',{obj.steady_state_blocks.variables},...
        'equations',{obj.steady_state_blocks.equations},...
        'swapping_func',obj.steady_state_blocks.swapping_func,...
        'extend_param_func',obj.steady_state_blocks.extend_param_func,...
        'is_log_var',obj.steady_state_blocks.is_log_var,...
        'solve_bgp',obj.options.solve_bgp,...
        'solve_bgp_shift',obj.options.solve_bgp_shift,...
        'sseqtns',{obj.routines.static.separate_blocks},...
        'ssjacobians',{obj.routines.static.separate_jacobians},...
        'debug',obj.options.debug,...
        'endo_list',{obj.endogenous.name},...
        'dynamic',{obj.equations.dynamic},...
        'residcode',obj.routines.static.one_block,...
        'auxcode',auxcode,...
        'arg_zero_solver',obj.options.steady_state_solver,...
        'optimopt',optimopt,...
        'solve_linear',obj.options.solve_linear);
    
    obj.steady_state_2_blocks_optimization=blocks;
    
else
    
    blocks=obj.steady_state_2_blocks_optimization;
    
end



end



function obj=recreate_steady_state_functions(obj)

% the steady state is always ordered alphabetically
%---------------------------------------------------
ss_occurrence=sum(obj.fast_sstate_occurrence.y,3);

ss_occurrence(ss_occurrence>0)=1;

swaps=obj.options.steady_state_endo_param_swap;

is_log_var=obj.endogenous.is_log_var;

if isempty(swaps)
    
    [swappingFunc,extendParamFunc]=swap_func();
    
else
    
    swaps_values=cell2mat(swaps(:,2));
    
    pos_endo=locate_variables(swaps(:,1),obj.endogenous.name);
    
    pos_param=locate_variables(swaps(:,3),obj.parameters.name);
    
    occparam=obj.fast_sstate_occurrence.param(:,pos_param);
    
    % swap occurrences: now the parameters take the place of the variables
    
    ss_occurrence(:,pos_endo)=occparam;
    
    % the following function re-establishes the order before each
    % evaluation
    
    [swappingFunc,extendParamFunc]=swap_func(pos_endo,pos_param,swaps_values);
    
    is_log_var(pos_endo)=false;
    
end

% nullify the columns of the variables that are fixed append the columns of
% the parameters that are endogenized...

if obj.options.solve_linear
    reblocks=false;
else
    reblocks=obj.options.solve_sstate_blocks;
end

[equations_blocks,variables_blocks]=parser.block_triangularize(ss_occurrence,...
    reblocks);

obj.steady_state_blocks=struct('equations',{equations_blocks},...
    'variables',{variables_blocks},...
    'swapping_func',swappingFunc,...
    'extend_param_func',extendParamFunc,...
    'is_log_var',is_log_var);

try
    static=obj.equations.shadow_fast_ssmodel;
catch
    static=obj.equations.shadow_static;
end

sssae=obj.equations.shadow_steady_state_auxiliary_eqtns;

repl_log=@replace_log; %#ok<NASGU>

repl_lin=@replace_lin; %#ok<NASGU>

if obj.options.solve_bgp
    
    is_log_var=obj.endogenous.is_log_var;
    
    endo_nbr=numel(is_log_var);
    
    % replace the time subscripts
    for ivar=1:endo_nbr
        
        digit= sprintf('%g',ivar);
        
        expr=['y\((',digit,')\)',...
            '\{(\+|\-)?(\d+)\}'];
        
        if is_log_var(ivar)
            
            static=regexprep(static,expr,'${repl_log($1,$2,$3)}');
            
            sssae=regexprep(sssae,expr,'${repl_log($1,$2,$3)}');
            
        else
            
            static=regexprep(static,expr,'${repl_lin($1,$2,$3)}');
            
            sssae=regexprep(sssae,expr,'${repl_lin($1,$2,$3)}');
            
        end
        
    end
    
    % auxiliary equation for solving for growth rates
    %-------------------------------------------------
    
    sssae_gr=auxil_growth_rates(sssae,is_log_var);
    
else
    
    % remove the time subscripts
    expr='(\{[\+\-]?\d+\})';
    
    static=regexprep(static,expr,'');
    
    sssae=regexprep(sssae,expr,'');
    
    sssae_gr={};
    
end

% create functions
list={'y','g','x','param','def'};

% TODO: Link this up with solve_function_mode: explicit, amateur,
% vectorized, professional
do_vectorize=false;

neqtns=numel(static);

ynames=utils.char.create_names([],'y',neqtns);

gnames=utils.char.create_names([],'g',neqtns);

static=utils.code.code2func(static,list,do_vectorize);

nblks=numel(equations_blocks);

theBlks=cell(1,nblks);

theJacobs=cell(1,nblks);

devectorize=true;

for iblk=1:nblks
    
    theBlks(iblk)=utils.code.code2vector(static(equations_blocks{iblk}),...
        devectorize);
    
    if obj.options.solve_linear||obj.options.steady_state_use_jacobian
        
        do_jacobian();
        
    end
    
end

OneBlk=utils.code.code2vector(static,devectorize);

obj.routines.static=struct('one_block',OneBlk{1},'separate_blocks',{theBlks},...
    'separate_jacobians',{theJacobs});

obj.routines.shadow_steady_state_auxiliary_eqtns=...
    struct('code',cell2mat(sssae(:)'),...
    'argins',{list},...
    'argouts',{{'y'}});

obj.routines.shadow_steady_state_auxiliary_growth_rate_eqtns=...
    struct('code',cell2mat(sssae_gr(:)'),...
    'argins',{list},...
    'argouts',{{'g'}});

    function []=do_jacobian()
        
        ywrt=ynames(variables_blocks{iblk});
        
        gwrt=[];
        
        order=1;
        
        if obj.options.solve_bgp
            
            gwrt=gnames(variables_blocks{iblk});
            
        end
        
        myfunc=static(equations_blocks{iblk});
        
        wrt=[ywrt(:).',gwrt(:).'];
        
        [derivs]=parser.differentiate_system(myfunc,list,wrt,order);
        
        %         if iblk==1
        %
        %             bad_fields=fieldnames(derivs)-derivative_fields;
        %
        %         end
        %
        %         theJacobs{iblk}=rmfield(derivs,bad_fields);
        
        theJacobs{iblk}=derivs;
        
    end


    function out=replace_lin(vp,sign_,digit)
        
        if isempty(sign_)||strcmp(sign_,'+')
            
            out=['(y(',vp,')+',digit,'*g(',vp,'))'];
            
        else
            
            out=['(y(',vp,')-',digit,'*g(',vp,'))'];
            
        end
        
    end


    function out=replace_log(vp,sign_,digit)
        
        if isempty(sign_)||strcmp(sign_,'+')
            
            out=['(y(',vp,')*g(',vp,')^',digit,')'];
            
        else
            
            out=['(y(',vp,')/g(',vp,')^',digit,')'];
            
        end
        
    end

end



function sssae=auxil_growth_rates(sssae,islogvar)

myreplace1=@just_neat1; %#ok<NASGU>

patt1='y\((\d+)\)=([^;]+);';

repl1='${myreplace1($1,$2)}';

sssae=regexprep(sssae,patt1,repl1);

    function out=just_neat1(vn,formule)
        
        myreplace2=@just_lag_it;  %#ok<NASGU>
        
        patt2='y\((\d+)\)';
        
        repl2='${myreplace2($1)}';
        
        newFormula=regexprep(formule,patt2,repl2);
        
        if islogvar(str2double(vn))
            % ratio
            out=['g(',vn,')=(',formule,')/(',newFormula,');'];
        else
            % subtract
            out=['g(',vn,')=',formule,'-(',newFormula,');'];
        end
        
    end

    function out=just_lag_it(vn)
        
        if islogvar(str2double(vn))
            % ratio
            out=['(y(',vn,')/g(',vn,'))'];
        else
            % subtract
            out=['(y(',vn,')-g(',vn,'))'];
        end
        
    end

end



function [sf,es]=swap_func(pos_endo,pos_param,swaps_values)

n=nargin;

sf=@swap_engine;

es=@extend_swaps;

    function [y,p]=swap_engine(y,p)
        
        if n>0
            
            p(pos_param)=y(pos_endo);
            
            y(pos_endo)=swaps_values;
            
        end
        
    end

    function p=extend_swaps(p,pold)
        
        if n>0
            
            ncols=size(p,2);
            
            p(pos_param,:)=pold(pos_param,ones(1,ncols));
            
        end
        
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
        
        % the steady state function may return both the steady state and
        % the balance growth path
        if any(isnan(varargout{1}(:)))||any(isinf(varargout{1}(:)))
            
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



function auxcode=auxiliary_memoizer(aux_ssfunc,aux_ss_gr_func,planner_routines,is_lagr_mult)

auxcode=@auxiliary_evaluation;

    function [yg,retcode]=auxiliary_evaluation(yg,x,p,d)
        
        y=yg(:,1);
        
        g=[];
        
        is_growth=size(yg,2)==2;
        
        if is_growth
            
            g=yg(:,2);
            
        end
        
        retcode=0;
        
        if ~isempty(aux_ssfunc) && (isa(aux_ssfunc,'function_handle')||...
                (isstruct(aux_ssfunc)&&~isempty(aux_ssfunc.code)))
            
            y=utils.code.evaluate_functions(aux_ssfunc,y,g,x,p,d);
            
            if is_growth
                
                g=utils.code.evaluate_functions(aux_ss_gr_func,y,g,x,p,d);
                
            end
            
        end
        
        if ~isempty(planner_routines)
            
            % TODO: Think of a way to compute the growth rate of the
            % lagrange multipliers besides doing it numerically through
            % divide and conquer.
            
            siz=planner_routines.planner_static_mult_support{1};
            
            pos=planner_routines.planner_static_mult_support{2};
            
            vals=zeros(size(pos));
            
            ss=y;
            
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
        
        yg=[y,g];
        
    end

end