function [Tz,eigval,retcode,options]=solve_first_order(structural_matrices,...
    siz,pos,options,k_future)
% INTERNAL FUNCTION
%

dv=structural_matrices.dv;

Q=structural_matrices.transition_matrices.Q;

[sm]=utils.solve.pull_first_order_partitions(dv,pos.v);

old_Tz=[];

if isfield(structural_matrices.old_solution,'Tz')
    
    old_Tz=structural_matrices.old_solution.Tz;
    
end

misalignment=structural_matrices.misalignment;

% the sm that comes out here is aggregated!
[Tz_pb,eigval,retcode,options,sm]=dsge_solver_first_order_autoregress_h(sm,Q,siz,pos,options,old_Tz);

nsols=size(Tz_pb,4);

Tz=cell(1,siz.h,nsols);

for isol=1:nsols
    
    one_solution(Tz_pb(:,:,:,isol))
    
end

    function one_solution(Tz_pb)
        
        if ~retcode
            % (non-)certainty equivalence
            %----------------------------
            dt_t=cell2mat(sm.sig_0); 
            
            if isempty(dt_t)
                
                dt_t=zeros(siz.nd,siz.h);
                
            end
            
            % add the residuals
            %-------------------
            dt_t=dt_t+structural_matrices.user_resids;
            
            xx=structural_matrices.y_ordered;
            
            gg=(1-options.solve_kill_g)*structural_matrices.bgp_ordered;
            
            A0sig=zeros(siz.nd,siz.nT,siz.h);
            
            Tz_sig=zeros(siz.nT,siz.h);
            
            df_Lf_Tzp_Lp=zeros(siz.nd,siz.nT);
            
            db_Lb_Tzb_Lb=zeros(siz.nd,siz.nT);
            
            UUi=zeros(siz.nd,siz.nT,siz.h);
                        
            Tz_e_rt=cell(1,siz.h);
            
            for rt=1:siz.h
                
                for rplus=1:siz.h
                    
                    % provision for shock impacts
                    %----------------------------
                    if ~options.occbin.do_it(rplus)
                        
                        continue
                        
                    end
                    
                    ds_0=dv{rt,rplus}(:,pos.v.s_0);
                    
                    dp_0=dv{rt,rplus}(:,pos.v.p_0);
                    
                    db_0=dv{rt,rplus}(:,pos.v.b_0);
                    
                    df_0=dv{rt,rplus}(:,pos.v.f_0);
                    
                    A0_0_1=[ds_0,dp_0,db_0,df_0];
                    
                    % provision for non-certainty equivalence
                    %----------------------------------------
                    df_plus=dv{rt,rplus}(:,pos.v.f_plus);
                    
                    db_plus=dv{rt,rplus}(:,pos.v.b_plus);
                    
                    dbf_plus=dv{rt,rplus}(:,pos.v.bf_plus);
                    
                    pb_=pos.t.pb;
                    
                    gTx=gg(:,rplus);
                    
                    if misalignment
                        
                        gTx=gTx+Tz_pb(:,:,rplus)*(xx(pb_,rplus)-xx(pb_,rt)-gg(pb_,rt));
                        
                    end
                    
                    dt_t(:,rt)=dt_t(:,rt)-dbf_plus*gTx(pos.t.bf,1);
                    
                    df_Lf_Tzp_Lp(:,pos.t.p)=df_plus*Tz_pb(pos.t.f,1:siz.np,rplus);% place in the p position
                    
                    db_Lb_Tzb_Lb(:,pos.t.b)=db_plus*Tz_pb(pos.t.b,siz.np+(1:siz.nb),rplus);% place in the b position
                    
                    A0sig(:,:,rt) = A0sig(:,:,rt) + A0_0_1 + df_Lf_Tzp_Lp + db_Lb_Tzb_Lb;
                    
                    % provision for shock impacts
                    %----------------------------
                    UUi(:,:,rt)=UUi(:,:,rt)+A0_0_1;
                    
                    UUi(:,pos.t.pb,rt)=UUi(:,pos.t.pb,rt)+sm.dbf_plus{rt,rplus}*Tz_pb(pos.t.bf,:,rplus);
                    
                end
                
                % shock impacts (current)
                %------------------------
                if ~options.occbin.do_it(rt)
                    
                    continue
                    
                end
                
                UUi(:,:,rt)=UUi(:,:,rt)\eye(siz.nT);
                
                Tz_e_rt{rt}=-UUi(:,:,rt)*sm.de_0{rt};
                
                Tz{1,rt,isol}=[Tz_pb(:,:,rt),Tz_sig(:,rt),Tz_e_rt{rt}];
                
            end
            
            if ~isempty(options.solve_occbin)
                
                options.occbin.B=sm.de_0;
                
                if k_future
                    
                    error('Occbin-type models not solved with anticipated events')
                    
                end
                
            end
            
            % shock impacts (future)
            %-----------------------
            for ik=1:k_future
                
                Tz_e_r0=Tz_e_rt;
                
                for rt=1:siz.h
                    
                    Sdbf_plus_rt_Tz_e0=0;
                    
                    for rplus=1:siz.h
                        
                        Sdbf_plus_rt_Tz_e0=Sdbf_plus_rt_Tz_e0+sm.dbf_plus{rt,rplus}*Tz_e_r0{rplus}(pos.t.bf,:);
                        
                    end
                    
                    Tz_e_rt{rt}=-UUi(:,:,rt)*Sdbf_plus_rt_Tz_e0;
                    
                    Tz{1,rt,isol}=[Tz{1,rt,isol},Tz_e_rt{rt}];
                    
                end
                
            end
            
            % now solve sum(A+*Tz_sig(+)+A0_sig*Tz_sig+dt_t=0
            %-------------------------------------------------
            % first we augment dt_t with the user_resid so that the first-order
            % approximation recoups the zero-th order approximation.
            
            Tz_sig=zeros(size(dt_t));
            
            if isempty(options.solve_occbin)
                
                if max(abs(dt_t(:)))>1e-9
                    
                    Tz_sig=dsge_tools.solve_perturbation_impact(Tz_sig,A0sig,...
                        sm.dbf_plus,dt_t,siz,pos);
                    
                end
                
            else
                
                options.occbin.c=dt_t;
                
            end
            
            if any(Tz_sig(:))
                
                for rt=1:siz.h
                    
                    Tz{1,rt,isol}(:,siz.np+siz.nb+1)=Tz_sig(:,rt);
                    
                end
                
            end
            % ensure the result is sparse
            %-----------------------------
            for rt=1:siz.h
                
                if ~isempty(Tz{1,rt,isol})
                    
                    Tz_proto=Tz{1,rt,isol};
                    
                    Tz{1,rt,isol}=sparse(Tz_proto);
                    
                end
                
            end
            
            for rt=1:siz.h
                
                if ~options.occbin.do_it(rt)
                    
                    Tz{1,rt,isol}=nan(size(Tz_proto));
                    
                end
                
            end
            
        end
        
    end

end