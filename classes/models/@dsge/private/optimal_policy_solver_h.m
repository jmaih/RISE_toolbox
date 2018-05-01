function [T,eigval,retcode,obj]=optimal_policy_solver_h(obj,structural_matrices,do_alphabet)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:


if isempty(obj)
    
    mydefaults=the_defaults();
    
    if nargout
        
        T=mydefaults;
        
    else
        
        disp_defaults(mydefaults);
        
    end
    
    return

end

if nargin<3
    
    do_alphabet=true;

end

gam=structural_matrices.planner.commitment{1};

beta=structural_matrices.planner.discount{1};

Q=structural_matrices.transition_matrices.Qinit;

shock_horizon=max(obj.exogenous.shock_horizon(:));

h=size(Q,1);

% form Aplus, A0, Aminus, B
%--------------------------
[sm]=utils.solve.pull_first_order_partitions(structural_matrices.dv,obj.locations.before_solve.v);

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
        
        a0_rt=a0_rt+[sm.ds_0{rt,rplus},sm.dp_0{rt,rplus},sm.db_0{rt,rplus},sm.df_0{rt,rplus}];
        
        aminus_rt=aminus_rt+sm.dpb_minus{rt,rplus};
        
        b_rt=b_rt+sm.de_0{rt,rplus};
        
        if rt==1 && rplus==1
            
            [neqtns,ny]=size(a0_rt);
            
            aplus=zeros(neqtns,ny);
            
            Aminus__=aplus;
        
        end
        
        aplus(:,bf_locs)=sm.dbf_plus{rt,rplus};
        
        Aplus{rt,rplus}=sparse(aplus);
    
    end
    
    A0{rt}=sparse(a0_rt);
    
    Aminus__(:,pb_locs)=aminus_rt;
    
    Aminus{rt}=sparse(Aminus__);
    
    B{rt}=sparse(b_rt);

end

clear aplus Aminus__

[AAplus,AA0,AAminus,BB,WW,reordering_index]=recreate_partial_system();

H0=[];

% solve the problem
%------------------
[H,G,retcode]=loose_commitment_engine(gam,beta,WW,AAplus,AA0,AAminus,BB,Q,...
    H0,shock_horizon,obj.options);

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
    
end

    function [AAplus,AA0,AAminus,BB,WW,reordering_index]=recreate_partial_system()
        
        if do_alphabet
            
            alphabetize()
        
        end
        
        [new_y_order,mult_cols,y_cols,orig_eqtns_rows,reordering_index]=load_help_indexes();
        % the order in the following call is strict and may not be changed
        %-----------------------------------------------------------------
        [AAplus,AA0,AAminus,BB,WW]=extract(Aplus,A0,Aminus,B);
        
        function [new_y_order,mult_cols,y_cols,orig_eqtns_rows,...
                reordering_index]=load_help_indexes()
            
            if ~isfield(obj.planner_system,'')
                % first put back in alphabetical order
                %-------------------------------------
                ylist=obj.endogenous.name;
                
                order_var_list=ylist(obj.order_var);
                
                if do_alphabet
                
                    mult_cols=obj.endogenous.is_lagrange_multiplier;
                
                else
                    
                    mult_cols=obj.endogenous.is_lagrange_multiplier(obj.order_var);
                
                    ylist=order_var_list;
                
                end
                
                y_cols=~mult_cols;
                
                mult_list=ylist(mult_cols);
                
                ylist=ylist(y_cols);
                
                % multipliers ordered last in the solution
                %-----------------------------------------
                
                grand_list=[ylist,mult_list];
                
                reordering_index=locate_variables(order_var_list,grand_list);
                
                new_y_order=locate_variables(ylist,obj.planner_system.wrt);
                
                nmult=sum(mult_cols);
                
                norigeqtns=nmult;
                
                % in the derivations, the original equations are placed on top
                %--------------------------------------------------------------
                
                orig_eqtns_rows=true(1,obj.endogenous.number);
                
                orig_eqtns_rows(norigeqtns+1:end)=false;
                
                obj.planner_system.new_y_order=new_y_order;
                
                obj.planner_system.mult_cols=mult_cols;
                
                obj.planner_system.y_cols=y_cols;
                
                obj.planner_system.orig_eqtns_rows=orig_eqtns_rows;
            
                obj.planner_system.reordering_index=reordering_index;
            
            end
            
            new_y_order=obj.planner_system.new_y_order;
            
            mult_cols=obj.planner_system.mult_cols;
            
            y_cols=obj.planner_system.y_cols;
            
            orig_eqtns_rows=obj.planner_system.orig_eqtns_rows;
            
            reordering_index=obj.planner_system.reordering_index;
        
        end

        function varargout=extract(varargin)
            
            nargs=numel(varargin);
            
            if nargout~=nargs+1
                
                error('mismatch between number of inputs and outputs')
            
            end
            
            varargout=cell(1,nargs+1);
            
            varargout(1:nargs)=varargin;
            
            for iarg=1:nargs
                
                arg=varargin{iarg};
                
                for ii=1:numel(arg)
                    
                    if iarg==nargs
                        
                        b=arg{ii}(orig_eqtns_rows,:);
                    
                    else
                        
                        b=arg{ii}(orig_eqtns_rows,y_cols);
                        
                        if obj.options.debug && iarg==2
                            % this works only for A0 since the coefficient
                            % matrices for the multipliers are swapped and
                            % scaled by powers of beta.
                            atest=-arg{ii}(~orig_eqtns_rows,mult_cols).';
                            
                            disp(max(max(abs(b-atest(:,new_y_order)))))
                            
                            keyboard
                        
                        end
                        
                    end
                    
                    varargout{iarg}{ii}=b;
                    % add the weights
                    %-----------------
                    if iarg==2
                        
                        varargout{end}{ii}=weights_extractor();
                    
                    end
                    
                end
                
            end
            
            function W=weights_extractor()
                % minus-ize so as to reproduce the same solution as the
                % non-loose commitment program.
                W=-.5*arg{ii}(~orig_eqtns_rows,y_cols);
                
                % re-order the rows according to the order in use
                W=W(new_y_order,:);
            
            end
            
        end
        
        function alphabetize()
            
            iov=obj.inv_order_var;
            
            A0=engine(A0);
            
            Aplus=engine(Aplus);
            
            Aminus=engine(Aminus);
            
            function b=engine(a)
                
                b=a;
                
                for ii=1:numel(a)
                    
                    b{ii}=a{ii}(:,iov);
                
                end
                
            end
            
        end
        
    end

    function re_order_solution()
        
        H=H(reordering_index,reordering_index,:);
        
        H=H(:,obj.locations.after_solve.t.pb,:);
        
        G=G(reordering_index,:,:,:);
    
    end

end

function [H,G,retcode]=loose_commitment_engine(gam,beta,W,Aplus,A0,Aminus,...
    B,Q,H0,k,options)

[nmult,ny]=size(A0{1});

h=size(Q,1);

nx=size(B{1},2);

n=nmult+ny;

y=1:ny;

lamb=ny+1:n;

lc_reconvexify=options.lc_reconvexify;

lc_use_pinv=options.lc_use_pinv;

GAMm=big_gam1();

if isempty(H0)
    
    switch options.solve_initialization
        
        case 'zeros'
            
            H0=zeros(n,n,h);
        
        case 'random'
            
            H0=randn(n,n,h);
        
        case 'backward'
            
            H0=backward_sol();
        
        otherwise
            
            error(['initialization ',parser.any2str(options.solve_initialization),' unknown'])
    
    end
    
    H0=iterate_func(H0);

end

H0=reshape(H0,n,n,h);

[H,~,retcode]=fix_point_iterator(@iterate_func,H0,options);

GAM0i=cell(1,h);

for rt=1:h
    
    if lc_use_pinv
        
        GAM0i{rt}=pinv(big_gam0(rt));
    
    else
        
        GAM0i{rt}=big_gam0(rt)\eye(n);
    
    end
    
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

    function [H1,H1_H0]=iterate_func(H00)
        
        H1=H00;
        
        H0=H00; % this is needed for the computation of GAM0
        
        for r0=1:h
            
            GAM0=big_gam0(r0);
            
            if lc_use_pinv
                
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
            
            if gam>0||lc_reconvexify
                
                tmp=0;
                % this is just an approximation
                
                for rplus=1:h
                    
                    tmp=tmp+Aplus{r0,rplus};
                
                end
                
                if lc_reconvexify
                    
                    GAM1(y,lamb)=gam/beta*tmp.';
                
                else
                    
                    GAM1(y,lamb)=1/beta*tmp.';
                
                end
                
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

    function H0=backward_sol()
        
        H0=zeros(n,n,h);
        
        AA0=zeros(n);
        
        AAminus=zeros(n);
        
        for rt_=1:h
            
            AA0(1:ny,1:ny)=2*W{rt_};
            
            AA0(1:ny,ny+1:end)=A0{rt_}.';
            
            AA0(ny+1:end,1:ny)=A0{rt_};
            
            AAminus(ny+1:end,1:ny)=Aminus{rt_};
            
            if lc_use_pinv
                
                H0(:,:,rt_)=pinv(AA0)*AAminus;
            
            else
                
                H0(:,:,rt_)=AA0\AAminus;
            
            end
            
        end
        
    end

end

function d=the_defaults()

d={
    'lc_algorithm','short',@(x)ismember(x,{'short','long'}),...
    'lc_algorithm must be a logical "short" or "long"'
    
        'lc_reconvexify',false,@(x)islogical(x),'lc_reconvexify must be a logical'
        
        'lc_use_pinv',true,@(x)islogical(x),'lc_use_pinv must be a logical'
        
        };
    
end