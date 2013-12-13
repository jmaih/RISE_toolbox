function [T,R,SS,retcode,optim_opt,itercode,algo_name]=dsge_lc_solve(...
    Aminus,A0,Aplus,B,W,gam,betta,...
    solve_expect_order,reordering_index,T0,optim_opt)

% Reference: Debortoli, Maih, Nunes (2010): "Loose Commitment in Medium
% and Large-Scale Models"
% This version, June 13, 2011

if nargin==0
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    T=loose_commitment_solver();
    return
end

if nargin<11
    optim_opt=[];
    if nargin<10
        T0=[];
        if nargin<9
            reordering_index=[];
            if nargin<8
                solve_expect_order=1;
                if nargin <7
                    betta=1;
                    if nargin<6
                        gam=1;
                        if nargin<5
                            error([mfilename,':: at least arguments Aminus,A0,Aplus,B,C,Q should be passed'])
                        end
                    end
                end
            end
        end
    end
elseif nargin>11
    error([mfilename,':: number of arguments cannot exceed 10'])
end

if isfield(optim_opt,'solve_accelerate')
    lc_accelerate_solver=optim_opt.solve_accelerate;
else
    lc_accelerate_solver=true;
end

[mult_nbr,endo_nbr]=size(A0);
m=endo_nbr+mult_nbr;
debug=false;
if lc_accelerate_solver
    exo_nbr=size(B{1},2);
    % find the static variables that are part of the objective function.
    enforce_no_static=find(any(W));
    [AA0,AAminus,AAplus,BB,stat_cols,QQ]=computational_savings(A0,Aminus,Aplus,enforce_no_static,B{:});
    orig_order=1:m;
    mult_cols=[false(1,endo_nbr),true(1,mult_nbr)];
    stat_cols=[stat_cols,false(1,mult_nbr)];
    dyn_cols=~stat_cols & ~mult_cols;
    
    static_endo_nbr=sum(stat_cols);
    dyamic_endo_nbr=sum(dyn_cols);
    
    if ~isempty(T0)
        % I need to reorder T0 and cut the part that goes into the reduced
        % system and I don't want to do that right now.
    end
    
    tmpBB=BB{1};
    for ireg=2:numel(BB)
        tmpBB(:,:,ireg)=BB{ireg};
    end
    clear BB
    [TT,RR,retcode,optim_opt,itercode,algo_name]=loose_commitment_solver(...
        AAminus(static_endo_nbr+1:end,static_endo_nbr+1:end),...
        AA0(static_endo_nbr+1:end,static_endo_nbr+1:end),...
        AAplus(static_endo_nbr+1:end,static_endo_nbr+1:end),tmpBB(static_endo_nbr+1:end,:,:),...
        W(dyn_cols(1:endo_nbr),dyn_cols(1:endo_nbr)),gam,betta,solve_expect_order,[],[],optim_opt);
    
    if retcode
        T=[]; R=[];
    else
        if any(stat_cols)
            % from the derivations, the order of the variables in the expanded T vector
            % will be xstat,xdyn,mu1,mu2. What we have got so far is T_xd_xd
            % T_xd_mu2bar,T_mu2_bar_xd,T_mu2_bar_mu2bar,R_xd, R_mu2bar
            
            x_r_d=1:dyamic_endo_nbr;
            mu2_r_d=dyamic_endo_nbr+(1:mult_nbr-static_endo_nbr);
            
            % new indices
            xs=1:static_endo_nbr;
            xd=static_endo_nbr+(1:dyamic_endo_nbr);
            mu1=static_endo_nbr+dyamic_endo_nbr+(1:static_endo_nbr);
            mu2=2*static_endo_nbr+dyamic_endo_nbr+1:m;
            %     x=[xs,xd];
            mu=[mu1,mu2];
            
            % extract the different parts
            Abar0_11=AA0(1:static_endo_nbr,1:static_endo_nbr);
            Abar0_12=AA0(1:static_endo_nbr,static_endo_nbr+1:end);
            Abar_minus_12=AAminus(1:static_endo_nbr,static_endo_nbr+1:end);
            Abar_plus_12=AAplus(1:static_endo_nbr,static_endo_nbr+1:end);
            Bbar_1=tmpBB(1:static_endo_nbr,:,:);
            clear AA0 AAminus AAplus 
            
            % initialize solution matrices
            nregs=size(tmpBB,3);
            T=zeros(m);
            R=zeros(m,exo_nbr,solve_expect_order,nregs);
            T(xd,xd)=TT(x_r_d,x_r_d);
            T(xd,mu2)=TT(x_r_d,mu2_r_d);
            T(mu2,xd)=TT(mu2_r_d,x_r_d);
            T(mu2,mu2)=TT(mu2_r_d,mu2_r_d);
            R(xd,:,:,:)=RR(x_r_d,:,:,:);
            R(mu2,:,:,:)=RR(mu2_r_d,:,:,:);
            clear  TT RR
            
            C=Abar_plus_12*T(xd,xd)+Abar0_12;
            iAbar0_11=Abar0_11\eye(static_endo_nbr);
            T(xs,xd)=-iAbar0_11*(C*T(xd,xd)+Abar_plus_12*T(xd,mu2)*T(mu2,xd)+Abar_minus_12);
            T(xs,mu2)=-iAbar0_11*(C*T(xd,mu2)+Abar_plus_12*T(xd,mu2)*T(mu2,mu2));
            
            for ireg=1:nregs
                R(xs,:,1,ireg)=-iAbar0_11*(C*R(xd,:,1,ireg)+Abar_plus_12*T(xd,mu2)*R(mu2,:,1,ireg)+Bbar_1(:,:,ireg));
                for ii=2:solve_expect_order
                    R(xs,:,ii,ireg)=-iAbar0_11*(C*R(xd,:,ii,ireg)+Abar_plus_12*T(xd,mu2)*R(mu2,:,ii,ireg)+Abar_plus_12*R(xd,:,ii-1,ireg));
                end
            end
            % the system has been solved using an auxiliary variable mu_tilde.
            % Before proceeding to solving for the static variables, we need to
            % express the system in terms of mu.
            T(xs,mu)=T(xs,mu)*QQ';
            T(xd,mu)=T(xd,mu)*QQ';
            T(mu,mu)=QQ*T(mu,mu)*QQ';
            T(mu,xd)=QQ*T(mu,xd);
            
            for ireg=1:nregs
                for h=1:solve_expect_order
                    R(mu,:,h,ireg)=QQ*R(mu,:,h,ireg);
                end
            end
            
            solve_order=[orig_order(stat_cols),orig_order(dyn_cols),orig_order(mult_cols)];
            T(solve_order,solve_order)=T;
            R(solve_order,:,:,:)=R;
        else
            T=TT;
            R=RR;
        end
        if debug
            % solve the problem without partialling out static variables and
            % compare
            [T_no_static,R_no_static]=loose_commitment_solver(...
                Aminus,A0,Aplus,B{:},W,gam,betta,solve_expect_order,[],[],optim_opt);
            keyboard
        end
    end
    
else
    tmpBB=B{1};
    for ireg=2:numel(B)
        tmpBB(:,:,ireg)=B{ireg};
    end
    [T,R,retcode,optim_opt,itercode,algo_name]=loose_commitment_solver(...
        Aminus,A0,Aplus,tmpBB,W,gam,betta,solve_expect_order,[],[],optim_opt);
end

if retcode
    SS=[];
else
    % finally re-order alphabetically
    if ~isempty(reordering_index)
        T=T(reordering_index,reordering_index);
        R=R(reordering_index,:,:,:);
    end
    SS=zeros(m,1);
% %     disp([mfilename,':: the steady state calculated here is not correct if there is a steady state file'])
end
end

