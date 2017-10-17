function [Tz_pb,eigval,retcode,options,sm]=...
    dsge_solver_first_order_autoregress_h(sm,Q,siz,pos,options,old_Tz)

% options
%--------
adjusted=struct();

adjusted.bf_cols=pos.t.bf;

adjusted.pb_cols=pos.t.pb;

adjusted.nd=siz.nd;

adjusted.siz=siz; % adjusted sizes

accelerate=options.solve_accelerate && siz.ns;

% aggregate A0 and A_
%--------------------
[sm,adjusted,accelSupport]=dsge_tools.aggregate_matrices(sm,siz,adjusted,accelerate);

if ~isempty(options.solve_occbin) && ~all(abs(diag(Q)-1)<1e-10)
    
    error('transition matrix must be diagonal for the occbin solution')
    
end

if isempty(options.solver)
    
    is_evs=is_eigenvalue_solver();
    
    if is_evs
        
        options.solver='rise_1';
        
    else
        
        options.solver='mfi';
        
    end
    
end

is_has_solution=~isempty(old_Tz);

if is_has_solution && accelerate
    % remove the solution for the static variables
    %----------------------------------------------
    for ireg=1:numel(old_Tz)
        
        old_Tz{ireg}=old_Tz{ireg}(siz.ns+1:end,:);
        
    end
    
end

T0=dsge_tools.utils.msre_initial_guess(sm.d0,sm.dpb_minus,sm.dbf_plus,...
    options.solve_initialization,old_Tz);

if is_has_solution
    
    Tz_pb=T0;
    
end

model_class=isempty(adjusted.pb_cols)+2*isempty(adjusted.bf_cols);

retcode=0;

eigval=[];

switch model_class
    
    case {1,3} % forward-looking or static models
        
        Tz_pb=if_then_else(is_has_solution,T0,0*T0);
        
    case 2 % backward-looking models
        
        Tz_pb=if_then_else(is_has_solution,T0,...
            dsge_tools.utils.msre_initial_guess(sm.d0,sm.dpb_minus,sm.dbf_plus,'backward',old_Tz));
        
    case 0 % hybrid models
        
        kron_method=strncmpi(options.solver,'mnk',3);
        
        is_known=true;
        
        if strcmpi(options.solver,'mfi')
            
            iterate_func=@(x)msre_solvers.functional_iteration_h(x,sm.dbf_plus,sm.d0,...
                sm.dpb_minus,adjusted.bf_cols,adjusted.pb_cols);
            
        elseif any(strcmpi(options.solver,{'mnk','mn'}))
            
            iterate_func=@(x)msre_solvers.newton_iteration_h(x,sm.dbf_plus,sm.d0,sm.dpb_minus,...
                adjusted.bf_cols,adjusted.pb_cols,kron_method,options);
            
        elseif any(strcmpi(options.solver,{'mfi_full','mnk_full','mn_full'}))
            
            [Gplus01,A0,Aminus,T0]=dsge_tools.full_state_matrices(adjusted.siz,sm,T0);
            
            if strcmpi(options.solver,'mfi_full')
                
                iterate_func=@(x)msre_solvers.functional_iteration_h_full(x,Gplus01,A0,Aminus);
                
            else
                
                iterate_func=@(x)msre_solvers.newton_iteration_h_full(x,Gplus01,A0,Aminus,...
                    kron_method,options);
                
            end
            
        elseif strcmpi(options.solver,'fwz')
            
            [Gplus01,A0,Aminus,T0]=dsge_tools.full_state_matrices(adjusted.siz,sm,T0);
            
            [iterate_func,solution_func,inverse_solution_func]= ...,sampling_func
                msre_solvers.fwz_newton_system(Gplus01,A0,Aminus,Q);
            
            T0=inverse_solution_func(T0);
            
        elseif any(strcmpi(options.solver,{'rise_1','klein','aim','sims'}))
            
            % don't do anything this is just so that does not flag the
            % solvers above
            
        else
            
            is_known=false;
            
        end
        
        if is_known
            
            switch lower(options.solver)
                
                case {'rise_1','klein','aim','sims'}
                    
                    dbf_plus_row=reconfigure_aplus();
                    
                    [Tz_pb,eigval,retcode]=dsge_solver_first_order_autoregress_1(...
                        dbf_plus_row,sm.ds_0,sm.dp_0,sm.db_0,sm.df_0,sm.dpb_minus,adjusted.siz,options);
                    
                    if ~isempty(options.solve_occbin)
                        
                        [options.occbin.Gplus01,options.occbin.A0,...
                            options.occbin.Aminus]=dsge_tools.full_state_matrices(adjusted.siz,sm,T0);
                        
                    end
                    
                case {'mfi','mfi_full','mnk','mnk_full','mn','mn_full','fwz'}
                    
                    [Tz_pb,~,retcode]=fix_point_iterator(iterate_func,T0,options);
                    
                    if  ~retcode && strcmpi(options.solver,'fwz')
                        
                        Tz_pb=solution_func(Tz_pb);
                        
                    end
                    
                    if any(strcmpi(options.solver,{'mfi_full','mnk_full','mn_full'}))
                        
                        Tz_pb=reshape(Tz_pb,[size(Tz_pb,1),size(Tz_pb,1),adjusted.siz.h]);
                        
                    end
                    
                    if any(strcmpi(options.solver,{'fwz','mfi_full','mnk_full','mn_full'}))
                        
                        Tz_pb=Tz_pb(:,adjusted.siz.ns+(1:adjusted.siz.np+adjusted.siz.nb),:);
                        
                    end
                    
            end
            
        else
            
            [user_solver,vargs]=utils.code.user_function_to_rise_function(...
                options.solver);
            
            % user-defined solver
            %--------------------
            [Gplus01,A0,Aminus,T0]=dsge_tools.full_state_matrices(adjusted.siz,sm,T0);
            
            [Tz_pb,~,retcode]=user_solver(Gplus01,A0,Aminus,Q,T0,...
                options.fix_point_TolFun,options.fix_point_maxiter,vargs{:});
            
            % collect the relevant part
            %--------------------------
            if ~retcode
                
                Tz_pb=Tz_pb(:,adjusted.siz.ns+(1:adjusted.siz.np+adjusted.siz.nb),:);
                
            end
            
        end
end


if ~retcode
    
    npb=siz.np+siz.nb;
    
    Tz_pb=reshape(Tz_pb,[adjusted.nd,npb,siz.h]);
    
    if accelerate
        % solve for the static variables
        %-------------------------------
        Sz_pb=zeros(siz.ns,npb,siz.h);
        
        for r0=1:siz.h
            
            ATT=0;
            
            for r1=1:siz.h
                
                ATT=ATT+accelSupport.Abar_plus_s{r0,r1}*Tz_pb(siz.np+1:end,:,r1); % <-- Tz_pb(npb+1:end,:,r1); we need also the both variables
                
            end
            
            ATT=ATT*Tz_pb(1:npb,:,r0);
            
            Sz_pb(:,:,r0)=-accelSupport.R_s_s{r0}\(accelSupport.Abar_minus_s{r0}+accelSupport.R_s_ns{r0}*Tz_pb(:,:,r0)+ATT);
            
        end
        
        Tz_pb=cat(1,Sz_pb,Tz_pb);
        
    end
    
end

    function Apl=reconfigure_aplus()
        
        Apl=cell(siz.h,1);
        
        for ii=1:siz.h
            
            if Q(ii,ii)
                
                Apl{ii}=sm.dbf_plus{ii,ii}/Q(ii,ii);
                
            else
                
                stud=find(Q(ii,:)>0,1,'first');
                
                Apl{ii}=sm.dbf_plus{ii,stud}/Q(ii,stud);
                % error('knife-edge probability matrix prevents from recovering a matrix')
            end
            
        end
        
    end

    function flag=is_eigenvalue_solver()
        
        flag=true;
        
        if ~(siz.h==1||all(diag(Q)==1))
            
            d0_test=sm.d0{1};
            
            dpb_minus_test=sm.dpb_minus{1};
            
            dbf_plus0=reconfigure_aplus();
            
            dbf_plus_test=dbf_plus0{1};
            
            % check whether it is shocks only
            for st=2:siz.h
                
                t0=get_max(d0_test-sm.d0{st});
                
                tplus=get_max(dbf_plus_test-dbf_plus0{st});
                
                tminus=get_max(dpb_minus_test-sm.dpb_minus{st});
                
                tmax=max([t0,tplus,tminus]);
                
                if tmax>1e-9
                    
                    flag=false;
                    
                    break
                    
                end
                
            end
            
        end
        
        function m=get_max(x)
            
            m=max(abs(x(:)));
            
        end
        
    end

end