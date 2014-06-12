function [T,eigval,zmultcols,retcode]=optimal_policy_solver_h(obj,structural_matrices)...

if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    T=struct('lc_reconvexify',false,...
        'lc_algorithm','short');
    return
end
gam=structural_matrices.planner.commitment{1};
beta=structural_matrices.planner.discount{1};
W=structural_matrices.planner.weights;
Q=structural_matrices.transition_matrices.Qinit;
shock_horizon=max(obj.exogenous.shock_horizon);
% re-order the rows of Aplus, A0, Aminus and B so that the lagrange
% multipliers are in the order : [static, pred, both, frwrd]
%-------------------------------------------------------------------
% eq_re_mult=equations_reordering_for_multipliers

reordering_index=obj.reordering_index;
eq_re_mult=obj.equations_reordering_for_multipliers;
h=size(Q,1);

% form Aplus, A0, Aminus, B
%--------------------------
[dbf_plus,ds_0,dp_0,db_0,df_0,dpb_minus,de_0]=...
    utils.solve.pull_first_order_partitions(structural_matrices.dv,...
    obj.locations.before_solve.v);

Aplus=cell(h);
A0=cell(1,h);
Aminus=cell(1,h);
B=cell(1,h);
bf_locs=obj.locations.before_solve.t.bf;
pb_locs=obj.locations.before_solve.t.pb;

for rt=1:h
    a0_rt=0;
    aminus_rt=0;
    b_rt=0;
    for rplus=1:h
        a0_rt=a0_rt+[ds_0{rt,rplus},dp_0{rt,rplus},db_0{rt,rplus},df_0{rt,rplus}];
        aminus_rt=aminus_rt+dpb_minus{rt,rplus};
        b_rt=b_rt+de_0{rt,rplus};
        if rt==1 && rplus==1
            [neqtns,ny]=size(a0_rt);
            aplus=zeros(neqtns,ny);
            Aminus__=aplus;
        end
        aplus(:,bf_locs)=dbf_plus{rt,rplus};
        Aplus{rt,rplus}=sparse(aplus);
    end
    A0{rt}=sparse(a0_rt);
    Aminus__(:,pb_locs)=aminus_rt;
    Aminus{rt}=sparse(Aminus__);
    B{rt}=sparse(b_rt);
end
clear aplus Aminus__

Aplus=re_order_rows(Aplus);
A0=re_order_rows(A0);
Aminus=re_order_rows(Aminus);
B=re_order_rows(B);
H0=[];

% solve the problem
%------------------
[H,G,zmultcols,retcode]=loose_commitment_engine(...
    gam,beta,W,Aplus,A0,Aminus,B,Q,H0,shock_horizon,obj.options);

T=[];
eigval=[];
if ~retcode
    % re-organize the solution in the order of order_var
    %---------------------------------------------------
    re_order_solution();
    
    % create final solution T
    %------------------------
    T=struct('Tz',{cell(1,h)});
    nrows=size(H,1);
    Tzsig=zeros(nrows,1);
    for ireg=1:h
        T.Tz{ireg}=[H(:,:,ireg),Tzsig,reshape(G(:,:,ireg,:),nrows,[])];
    end
    % augment multcols accordingly
    %-----------------------------
    exo_nbr=size(G,2);
    zmultcols=[zmultcols,false(1,exo_nbr*(shock_horizon+1))];
end

    function re_order_solution()
        H=H(reordering_index,reordering_index,:);
        H=H(:,obj.locations.after_solve.t.pb,:);
        G=G(reordering_index,:,:,:);
        zmultcols=zmultcols(reordering_index);
        % keep only the state multipliers
        %--------------------------------
        zmultcols=zmultcols(obj.locations.after_solve.t.pb);
    end
    function X=re_order_rows(X)
        for ii=1:numel(X)
            X{ii}=X{ii}(eq_re_mult,:);
        end
    end
end

function [H,G,multcols,retcode]=loose_commitment_engine(...
    gam,beta,W,Aplus,A0,Aminus,B,Q,H0,k,options)

[nmult,ny]=size(A0{1});
h=size(Q,1);
nx=size(B{1},2);
n=nmult+ny;
y=1:ny;
lamb=ny+1:n;
multcols=[false(1,ny),true(1,nmult)];

GAMm=big_gam1();

if isempty(H0)
    H0=zeros(n,n,h);
    use_pinv=true;
    H0=iterate_func(H0,use_pinv);
end

H0=reshape(H0,n,n,h);

[H,~,retcode]=fix_point_iterator(@iterate_func,H0,options);

GAM0i=cell(1,h);
for rt=1:h
    GAM0i{rt}=big_gam0(rt)\eye(n);
end

G=zeros(n,nx,h,k+1);
for ik=1:k+1
    if ik==1
        GAMe=big_gam_e0();
    else
        GAMe=big_gam_e_update(G(:,:,:,ik-1));
    end
    G(:,:,:,ik)=solve_shock_impact(GAMe);
end

    function [H1,H1_H0]=iterate_func(H00,use_pinv)
        if nargin<2
            use_pinv=false;
        end
        H1=H00;
        H0=H00; % this is needed for the computation of GAM0
        for r0=1:h
            GAM0=big_gam0(r0);
            if use_pinv
                H1(:,:,r0)=-pinv(GAM0)*GAMm{r0};
            else
                H1(:,:,r0)=-GAM0\GAMm{r0};
            end
        end
        H1_H0=H1-H00;
    end

    function Ge=solve_shock_impact(GAMe)
        Ge=G(:,:,:,1);
        for r0=1:h
            Ge(:,:,r0,1)=-GAM0i{r0}*GAMe(:,:,r0);
        end
    end

    function GAMe=big_gam_e_update(Ge_prev)
        GAMe=zeros(n,nx,h);
        for r0=1:h
            tmp=0;
            for rplus=1:h
                if gam>0
                    tmp=tmp+Q(r0,rplus)*Aminus{rplus}.'*Ge_prev(lamb,:,rplus);
                end
                GAMe(lamb,:,r0)=GAMe(lamb,:,r0)+Aplus{r0,rplus}*Ge_prev(y,:,rplus);
            end
            if gam>0
                GAMe(y,:,r0)=beta*gam*tmp;
            end
        end
    end

    function GAMe=big_gam_e0()
        GAMe=zeros(n,nx,h);
        for r0=1:h
            GAMe(lamb,:,r0)=B{r0};
        end
    end

    function GAMm=big_gam1()
        GAM1=zeros(n);
        GAMm=cell(1,h);
        for r0=1:h
            GAM1(lamb,y)=Aminus{r0};
            if gam>0
                tmp=0;
                % this is just an approximation
                for rplus=1:h
                    tmp=tmp+Aplus{r0,rplus};
                end
                GAM1(y,lamb)=1/beta*tmp.';
            end
            GAMm{r0}=GAM1;
        end
    end

    function GAM0=big_gam0(rt)
        GAM0=zeros(n);
        g0_y_l_1=0;
        g0_y_l_2=0;
        for rplus=1:h
            GAM0(y,y)=GAM0(y,y)+Q(rt,rplus)*Aminus{rplus}'*H0(lamb,y,rplus);
            
            g0_y_l_1=g0_y_l_1+Aplus{rt,rplus}*H0(y,y,rplus);
            if gam>0
                g0_y_l_2=g0_y_l_2+Q(rt,rplus)*Aminus{rplus}'*H0(lamb,lamb,rplus);
                GAM0(lamb,lamb)=GAM0(lamb,lamb)+Aplus{rt,rplus}*H0(y,lamb,rplus);
            end
        end
        GAM0(y,y)=2*W{rt}+beta*GAM0(y,y);
        
        GAM0(y,lamb)=(A0{rt}+(1-gam)*g0_y_l_1).';
        
        if gam>0
            GAM0(y,lamb)=GAM0(y,lamb)+beta*gam*g0_y_l_2;
            GAM0(lamb,lamb)=gam*GAM0(lamb,lamb);
        end
        
        GAM0(lamb,y)=A0{rt}+g0_y_l_1;
    end
end