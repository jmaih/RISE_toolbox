function [Tz_pb,eigval,retcode,options]=...
    dsge_solver_first_order_autoregress_h(dbf_plus,ds_0,dp_0,db_0,df_0,...
    dpb_minus,Q,siz,pos,options)

% options
%--------
bf_cols_adjusted=pos.t.bf;
pb_cols_adjusted=pos.t.pb;

nd_adjusted=siz.nd;
siz_adj=siz; % adjusted sizes
accelerate=options.solve_accelerate && siz.ns;
% aggregate A0 and A_
%--------------------
d0=cell(1,siz.h);
for r0=1:siz.h
    for r1=2:siz.h
        ds_0{r0,1}=ds_0{r0,1}+ds_0{r0,r1};
        dp_0{r0,1}=dp_0{r0,1}+dp_0{r0,r1};
        db_0{r0,1}=db_0{r0,1}+db_0{r0,r1};
        df_0{r0,1}=df_0{r0,1}+df_0{r0,r1};
        dpb_minus{r0,1}=dpb_minus{r0,1}+dpb_minus{r0,r1};
    end
    d0{r0}=[ds_0{r0,1},dp_0{r0,1},db_0{r0,1},df_0{r0,1}];
    % eliminate static variables for speed
    %-------------------------------------
    if accelerate
        if r0==1
            Abar_minus_s=cell(1,siz.h);
            R_s_s=cell(1,siz.h);
            R_s_ns=cell(1,siz.h);
            Abar_plus_s=cell(siz.h);
            nd_adjusted=nd_adjusted-siz.ns;
            bf_cols_adjusted=bf_cols_adjusted-siz.ns;
            pb_cols_adjusted=pb_cols_adjusted-siz.ns;
            siz_adj.ns=0;
            siz_adj.nd=siz_adj.nd-siz.ns;
            siz_adj.nT=siz_adj.nT-siz.ns;
        end
        [Q0,d0{r0}]=qr(d0{r0});
        dpb_minus{r0,1}=Q0'*dpb_minus{r0,1};
        Abar_minus_s{r0}=dpb_minus{r0,1}(1:siz.ns,:);
        dpb_minus{r0,1}=dpb_minus{r0,1}(siz.ns+1:end,:);
        for r1=1:siz.h
            dbf_plus{r0,r1}=Q0'*dbf_plus{r0,r1};
            Abar_plus_s{r0,r1}=dbf_plus{r0,r1}(1:siz.ns,:);
            dbf_plus{r0,r1}=dbf_plus{r0,r1}(siz.ns+1:end,:);
        end
        R_s_s{r0}=d0{r0}(1:siz.ns,1:siz.ns);
        R_s_ns{r0}=d0{r0}(1:siz.ns,siz.ns+1:end);
        d0{r0}=d0{r0}(siz.ns+1:end,siz.ns+1:end);
        ds_0{r0,1}=d0{r0}(:,1:siz_adj.ns);
        dp_0{r0,1}=d0{r0}(:,siz_adj.ns+(1:siz_adj.np));
        db_0{r0,1}=d0{r0}(:,siz_adj.ns+siz_adj.np+(1:siz_adj.nb));
        df_0{r0,1}=d0{r0}(:,siz_adj.ns+siz_adj.np+siz_adj.nb+(1:siz_adj.nf));
    end
end
ds_0=ds_0(:,1)';
dp_0=dp_0(:,1)';
db_0=db_0(:,1)';
df_0=df_0(:,1)';

if options.solve_occbin && ~all(abs(diag(Q)-1)<1e-10)
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

T0=dsge_tools.utils.msre_initial_guess(d0,dpb_minus,dbf_plus,...
    options.solve_initialization);

model_class=isempty(pb_cols_adjusted)+2*isempty(bf_cols_adjusted);

retcode=0;
eigval=[];
switch model_class
    case 1 % forward-looking models
        Tz_pb=0*T0;
    case 2 % backward-looking models
        Tz_pb=dsge_tools.utils.msre_initial_guess(d0,dpb_minus,dbf_plus,'backward');
    case 3 % static models
        Tz_pb=0*T0;
    case 0 % hybrid models
        kron_method=strncmpi(options.solver,'mnk',3);
        
        if strcmpi(options.solver,'mfi')
            iterate_func=@(x)msre_solvers.functional_iteration_h(x,dbf_plus,d0,...
                dpb_minus,bf_cols_adjusted,pb_cols_adjusted);
        elseif any(strcmpi(options.solver,{'mnk','mn'}))
            iterate_func=@(x)msre_solvers.newton_iteration_h(x,dbf_plus,d0,dpb_minus,...
                bf_cols_adjusted,pb_cols_adjusted,kron_method,options);
        elseif any(strcmpi(options.solver,{'mfi_full','mnk_full','mn_full'}))
            [Gplus01,A0,Aminus,T0]=full_state_matrices(siz_adj,dbf_plus,d0,dpb_minus,T0);
            if strcmpi(options.solver,'mfi_full')
                iterate_func=@(x)msre_solvers.functional_iteration_h_full(x,Gplus01,A0,Aminus);
            else
                iterate_func=@(x)msre_solvers.newton_iteration_h_full(x,Gplus01,A0,Aminus,...
                    kron_method,options);
            end
        elseif strcmpi(options.solver,'fwz')
            [Gplus01,A0,Aminus,T0]=full_state_matrices(siz_adj,dbf_plus,d0,dpb_minus,T0);
            [iterate_func,solution_func,inverse_solution_func]= ...,sampling_func
                msre_solvers.fwz_newton_system(Gplus01,A0,Aminus,Q);
            T0=inverse_solution_func(T0);
        end
        
        switch lower(options.solver)
            case {'rise_1','klein','aim','sims'}
                dbf_plus_row=reconfigure_aplus();
                [Tz_pb,eigval,retcode]=dsge_solver_first_order_autoregress_1(...
                    dbf_plus_row,ds_0,dp_0,db_0,df_0,dpb_minus,siz_adj,options);
                if options.solve_occbin
                    [Gplus01,A0,Aminus]=full_state_matrices(siz_adj,dbf_plus,d0,dpb_minus,T0);
                    options.occbin=struct('Gplus01',Gplus01,'A0',A0,'Aminus',Aminus);
                end
            case {'mfi','mfi_full','mnk','mnk_full','mn','mn_full','fwz'}
                [Tz_pb,~,retcode]=fix_point_iterator(iterate_func,T0,options);
                if  ~retcode && strcmpi(options.solver,'fwz')
                    Tz_pb=solution_func(Tz_pb);
                end
                if any(strcmpi(options.solver,{'mfi_full','mnk_full','mn_full'}))
                    Tz_pb=reshape(Tz_pb,[size(Tz_pb,1),size(Tz_pb,1),siz_adj.h]);
                end
                if any(strcmpi(options.solver,{'fwz','mfi_full','mnk_full','mn_full'}))
                    Tz_pb=Tz_pb(:,siz_adj.ns+(1:siz_adj.np+siz_adj.nb),:);
                end
            otherwise
                % user-defined solver
                %--------------------
                [Gplus01,A0,Aminus,T0]=full_state_matrices(siz_adj,dbf_plus,d0,dpb_minus,T0);
                
                [Tz_pb,~,retcode]=options.solver(Gplus01,A0,Aminus,Q,T0);
                
                % collect the relevant part
                %--------------------------
                if ~retcode
                    Tz_pb=Tz_pb(:,siz_adj.ns+(1:siz_adj.np+siz_adj.nb),:);
                end
        end
end


if ~retcode
    npb=siz.np+siz.nb;
    Tz_pb=reshape(Tz_pb,[nd_adjusted,npb,siz.h]);
    if accelerate
        % solve for the static variables
        %-------------------------------
        Sz_pb=zeros(siz.ns,npb,siz.h);
        for r0=1:siz.h
            ATT=0;
            for r1=1:siz.h
                ATT=ATT+Abar_plus_s{r0,r1}*Tz_pb(siz.np+1:end,:,r1); % <-- Tz_pb(npb+1:end,:,r1); we need also the both variables
            end
            ATT=ATT*Tz_pb(1:npb,:,r0);
            Sz_pb(:,:,r0)=-R_s_s{r0}\(Abar_minus_s{r0}+R_s_ns{r0}*Tz_pb(:,:,r0)+ATT);
        end
        Tz_pb=cat(1,Sz_pb,Tz_pb);
    end
end
    function Apl=reconfigure_aplus()
        Apl=cell(siz.h,1);
        for ii=1:siz.h
            if Q(ii,ii)
                Apl{ii}=dbf_plus{ii,ii}/Q(ii,ii);
            else
                error('knife-edge probability matrix prevents from recovering a matrix')
            end
        end
    end

    function flag=is_eigenvalue_solver()
        flag=true;
        if ~(siz.h==1||all(diag(Q)==1))
            d0_test=d0{1};
            dpb_minus_test=dpb_minus{1};
            dbf_plus0=reconfigure_aplus();
            dbf_plus_test=dbf_plus0{1};
            % check whether it is shocks only
            for st=2:siz.h
                t0=get_max(d0_test-d0{st});
                tplus=get_max(dbf_plus_test-dbf_plus0{st});
                tminus=get_max(dpb_minus_test-dpb_minus{st});
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