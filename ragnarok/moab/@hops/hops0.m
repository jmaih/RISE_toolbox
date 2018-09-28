classdef hops
    
    properties
        post
        posv
        posz
        n
        k
        lambda
        m
        h
        sig_in_v
        max_order=5
    end
    
    properties
        EU_Ee_links
        EU
        fv
        fvv
        fvvv
        fvvvv
        fvvvvv
        fplus_lambda_bf
    end
    
    properties
        ss
        g
        Tz
        Tzz
        Tzzz
        Tzzzz
        Tzzzzz
    end
    
    properties
        a0z
        a1z
        a0zz
        a1zz
        a0zzz
        a1zzz
        a0zzzz
        a1zzzz
        Jz
        Jzz
        Jzzz
        Jzzzz
        Jzzzzz
        B
    end
    
    properties
        kronall=@utils.kronecker.kronall
        sp=@utils.kronecker.sum_permutations
        cr2=@utils.cr.second_order
        cr3=@utils.cr.third_order
        cr4=@utils.cr.fourth_order
        cr5=@utils.cr.fifth_order
        keep
        expand
        Comp_mat
        UnComp_mat
        gensylv_solve=1
        debug
        use_kron=false
    end
    
    methods
        
        function obj=hops(ns,np,nb,nf,ne0,k,h,sig_in_v,debug_)
            
            if nargin
                
                if nargin<9
                    
                    debug_=false;
                    
                end
                
                obj.debug=debug_;
                
                obj.sig_in_v=sig_in_v;
                
                obj.k=k;
                
                obj.h=h;
                
                obj.n=set_sizes();
                
                [obj.post,obj.lambda]=tpositions();
                
                obj.posv=vpositions();
                
                [obj.posz,obj.m]=zpositions();
                
                [obj.keep,obj.expand,obj.Comp_mat,obj.UnComp_mat]=initialize_compress_expand();
                
                %------------------
                ch=cell(h);
                obj.a0z=ch;
                obj.a1z=ch;
                obj.a0zz=ch;
                obj.a1zz=ch;
                obj.a0zzz=ch;
                obj.a1zzz=ch;
                obj.a0zzzz=ch;
                obj.a1zzzz=ch;
                obj.Jz=ch;
                obj.Jzz=ch;
                obj.Jzzz=ch;
                obj.Jzzzz=ch;
                obj.Jzzzzz=ch;
                %--------------------
                
                obj.EU_Ee_links=set_EU_Ee_links(obj.n.e,obj.n.z,obj.m.sig,...
                    obj.max_order,obj.debug);
                
                [obj.EU,obj.max_order]=set_stochastic_matrices_moments(obj.n.e,...
                    obj.n.z,obj.max_order,obj.EU_Ee_links);
                
            end
            
            function [post,lambda]=tpositions()
                
                post=struct();
                
                [post.s,offset]=set_position(ns);
                
                [post.p,offset]=set_position(np,offset);
                
                [post.b,offset]=set_position(nb,offset);
                
                post.f=set_position(nf,offset);
                
                post.pb=[post.p,post.b];
                
                post.bf=[post.b,post.f];
                
                tmp=speye(obj.n.t);
                
                lambda=struct();
                
                lambda.pb=tmp(post.pb,:);
                
                lambda.bf=tmp(post.bf,:);
                
            end
            
            function posv=vpositions()
                
                posv=struct();
                
                [posv.bf,offset]=set_position(numel(obj.post.bf));
                
                [posv.spbf,offset]=set_position(obj.n.t,offset);
                
                [posv.pb,offset]=set_position(numel(obj.post.pb),offset);
                
                [posv.eps,offset]=set_position(obj.n.e,offset);
                
                [posv.sig,nv]=set_position(sig_in_v,offset);
                
                if nv~=obj.n.v
                    
                    error('mismatch in sizes')
                    
                end
                
            end
            
            function [posz,m]=zpositions()
                
                posz=struct();
                
                [posz.pb,offset]=set_position(numel(obj.post.pb));
                
                [posz.sig,offset]=set_position(1,offset);
                
                [posz.e0,offset]=set_position(obj.n.e,offset);
                
                [posz.e1_ek,nz]=set_position(obj.n.e*obj.k,offset);
                
                if nz~=obj.n.z
                    
                    error('mismatch in sizes')
                    
                end
                
                m=struct();
                
                tmp=speye(nz);
                
                m.pb=tmp(posz.pb,:);
                
                m.sig=tmp(posz.sig,:);
                
                m.e0=tmp(posz.e0,:);
                
                m.e1_ek=tmp(posz.e1_ek,:);
                
                m.e0_ek=[m.e0;m.e1_ek];
                
            end
            
            function [p,offset]=set_position(n,offset)
                
                % n can be logical and so it has to be doubled before
                % evaluation
                n=double(n);
                
                if nargin<2
                    
                    offset=0;
                    
                end
                
                p=offset+(1:n);
                
                if ~isempty(p)
                    
                    offset=p(end);
                    
                end
                
            end
            
            function n=set_sizes()
                
                n=struct();
                
                n.s=ns;
                
                n.p=np;
                
                n.b=nb;
                
                n.f=nf;
                
                n.e=ne0;
                
                n.t=ns+np+nb+nf;
                
                n.pb=np+nb;
                
                n.bf=nb+nf;
                
                n.z=n.pb+1+n.e*(obj.k+1);
                
                n.v=n.bf+n.t+n.pb+n.e+sig_in_v;
                
            end
            
        end
        
        function obj=solve(obj,f,f_iscompressed,Tz,ss,g,order,M)
            
            if nargin<8
                
                M=[];
                
            end
            
            absorb_derivatives()
            
            obj.Tz=Tz;
            
            obj.ss=ss;
            
            obj.g=g;
            
            if ~isempty(M)
                % reset the moments if the user specifies other moments
                % including e.g. skewed shocks
                obj.EU=set_stochastic_matrices_moments(obj.n.e,obj.n.z,order,...
                    obj.EU_Ee_links,M);
                
            end
            
            if order>1
                
                if obj.debug,tic,end
                
                obj=solve_second_order(obj);
                
                if obj.debug
                    
                    fprintf(1,'order 2 done in %0.4f\n',toc);
                    
                end
                
                if order>2
                    
                    if obj.debug,tic,end
                    
                    obj=solve_third_order(obj);
                    
                    if obj.debug
                        
                        fprintf(1,'order 3 done in %0.4f\n',toc);
                        
                    end
                    
                    if order>3
                        
                        if obj.debug,tic,end
                        
                        obj=solve_fourth_order(obj);
                        
                        if obj.debug
                            
                            fprintf(1,'order 4 done in %0.4f\n',toc);
                            
                        end
                        
                        if order>4
                            
                            if obj.debug,tic,end
                            
                            obj=solve_fifth_order(obj);
                            
                            if obj.debug
                                
                                fprintf(1,'order 5 done in %0.4f\n',toc);
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
            function absorb_derivatives()
                
                v=repmat('v',1,obj.max_order);
                
                true_order=min(obj.max_order,order);
                    
                obj=set_compress_uncompress(obj,true_order);
                
                for ii=1:true_order
                    
                    vi=v(1:ii);
                    
                    if ~isfield(f,vi)
                        
                        error(['derivatives of order ',int2str(ii),' not found'])
                        
                    end
                    
                    mykeep=obj.keep.v{ii};
                    
                    obj.(['f',vi])=store_derivatives(f.(vi),mykeep,f_iscompressed);
                    
                end
                
                clear f
                
            end
            
        end
        
        function obj=solve_second_order(obj)
            
            pb=obj.post.pb;
            
            bf=obj.post.bf;
            
            msig=obj.m.sig;
            
            sigpos=find(msig);
            
            mpb=obj.m.pb;
            
            me0=obj.m.e0;
            
            Jz_bulk=[msig
                obj.m.e1_ek
                sparse(obj.n.e,obj.n.z)
                ];
            
            a1z_bulk=sparse(obj.n.t+obj.n.pb+obj.n.e+obj.sig_in_v,obj.n.z);
            
            x=obj.ss;
            
            g_=obj.g;
            
            A=repmat({sparse(0)},1,obj.h);
            
            B_=A;
            
            CBAR=cell(obj.h);
            
            K=obj.Comp_mat.z{2};
            
            W=obj.UnComp_mat.z{2};
            
            for rt=1:obj.h
                
                for rtp1=1:obj.h
                    
                    obj.Jz{rt,rtp1}=do_Jz();
                    
                    %--------------------------------
                    Jz_=obj.Jz{rt,rtp1};
                    
                    %--------------------------------
                    obj.a0z{rt,rtp1}=do_a0z();
                    
                    obj.a1z{rt,rtp1}=do_a1z();
                    %--------------------------------                    
                    % #
                    a0z_=obj.a0z{rt,rtp1};
                    % #
                    a1z_=obj.a1z{rt,rtp1};
                    %--------------------------------
                    
                    A{rt}=A{rt}+a_iterate();
                    
                    B_{rt}=B_{rt}+b_iterate();
                    
                    % Compressed C
                    CBAR{rt,rtp1}=c_iterate();
                    
                end
                % compress
                A{rt}=A{rt}*K;
                
            end
            
            obj.B=B_; clear B_
            
            obj.fplus_lambda_bf=fplus_times_lambda_bf(obj.fv,obj.posv.bf,...
                obj.lambda.bf);
            
            obj.Tzz=solve_generalized_sylvester(A,obj.B,CBAR,...
                obj.fplus_lambda_bf,W,obj.gensylv_solve);
            
            function ai=a_iterate()
                
                fvv_=obj.fvv{rt,rtp1}(:,obj.expand.v{2});
                
                ai=fvv_Evz_vz();
                
                function r=fvv_Evz_vz()
                    
                    r=A_times_kron_Q1_Qk(obj.debug,fvv_,a0z_,a0z_);
                    
                    r=r+A_times_kron_Q1_Qk(obj.debug,fvv_,a1z_,a1z_)*obj.EU{2};
                    
                end
                
            end
            
            function bi=b_iterate()
                
                fplus=obj.fv{rt,rtp1}(:,obj.posv.bf);
                
                f0=obj.fv{rt,rtp1}(:,obj.posv.spbf);
                
                Tz1=obj.Tz{rtp1}(obj.post.bf,obj.posz.pb)*obj.lambda.pb;
                
                bi=fplus*Tz1+f0;
                
            end
            
            function ci=c_iterate()
                
                ci=A_times_kron_Q1_Qk(obj.debug,W,Jz_,Jz_)*K+...
                    W*obj.EU{2}*K;
                                
            end
            
            function Jz=do_Jz()
                
                x01=x(:,rtp1)-x(:,rt)-g_(:,rt);
                
                tmp=obj.Tz{rt};
                
                tmp(:,sigpos)=tmp(:,sigpos)-x01;
                
                Jz=[
                    tmp(pb,:)
                    Jz_bulk
                    ];
            end
            
            function a0z=do_a0z()
                
                TJg=obj.Tz{rtp1}*obj.Jz{rt,rtp1};
                
                TJg(:,sigpos)=TJg(:,sigpos)-g_(:,rtp1);
                
                if obj.sig_in_v
                    
                    a0z=[
                        TJg(bf,:)
                        obj.Tz{rt}
                        mpb
                        me0
                        msig
                        ];
                    
                else
                    
                    a0z=[
                        TJg(bf,:)
                        obj.Tz{rt}
                        mpb
                        me0
                        ];
                    
                end
                
            end
            
            function a1z=do_a1z()
                
                a1z=[obj.Tz{rtp1}(bf,:)
                    a1z_bulk];
                
            end
            
            function fv=fplus_times_lambda_bf(fv,posv_bf,lambda_bf)
                
                for ii=1:numel(fv)
                    
                    fv{ii}=fv{ii}(:,posv_bf)*lambda_bf;
                    
                end
                
            end

        end
        
        function obj=solve_third_order(obj)
            
            pb=obj.post.pb;
            
            bf=obj.post.bf;
            
            Jzz_bulk=sparse(obj.n.e*(obj.k+1)+1,obj.n.z^2);
            
            a1zz_bulk=sparse(obj.n.t+obj.n.pb+obj.n.e+obj.sig_in_v,obj.n.z^2);
            
            a0zz_bulk=sparse(obj.n.pb+obj.n.e+obj.sig_in_v,obj.n.z^2);
            
            omg1=utils.cr.omega(obj.n.z,1);
            
            A=repmat({sparse(0)},1,obj.h);
            
            CBAR=cell(obj.h);
            
            K=obj.Comp_mat.z{3};
            
            W=obj.UnComp_mat.z{3};
            
            for rt=1:obj.h
                
                for rtp1=1:obj.h
                    
                    obj.Jzz{rt,rtp1}=do_Jzz();
                    
                    %--------------------------------
                    Jz_=obj.Jz{rt,rtp1};
                    
                    Jzz_=obj.Jzz{rt,rtp1};
                    %--------------------------------
                    
                    obj.a0zz{rt,rtp1}=do_a0zz();
                    
                    obj.a1zz{rt,rtp1}=do_a1zz();
                    %--------------------------------                    
                    % #
                    a0z_=obj.a0z{rt,rtp1};
                    
                    a0zz_=obj.a0zz{rt,rtp1};
                    % #
                    a1z_=obj.a1z{rt,rtp1};
                    
                    a1zz_=obj.a1zz{rt,rtp1};
                    %--------------------------------
                    
                    A{rt}=A{rt}+a_iterate();
                    
                    CBAR{rt,rtp1}=c_iterate();
                    
                end
                % compress
                A{rt}=A{rt}*K;
                
            end
            
            obj.Tzzz=solve_generalized_sylvester(A,obj.B,CBAR,...
                obj.fplus_lambda_bf,W,obj.gensylv_solve);
            
            function Jzz=do_Jzz()
                
                if rtp1==1
                    
                    tmp=obj.Tzz{rt};
                    
                    Jzz=[
                        tmp(pb,:)
                        Jzz_bulk
                        ];
                    
                else
                    
                    Jzz=obj.Jzz{rt,1};
                    
                end
                
            end
            
            function a0zz=do_a0zz()
                
                Tzz_=obj.Tzz{rtp1};
                
                Tz_=obj.Tz{rtp1};
                
                TJg=obj.cr2(Tzz_,Tz_,Jzz_,Jz_);
                
                a0zz=[
                    TJg(bf,:)
                    obj.Tzz{rt}
                    a0zz_bulk
                    ];
            end
            
            function a1zz=do_a1zz()
                
                a1zz=[obj.Tzz{rtp1}(bf,:)
                    a1zz_bulk];
                
            end
            
            function ai=a_iterate()
                
                fplus_lmbd=obj.fv{rt,rtp1}(:,obj.posv.bf)*obj.lambda.bf;
                
                ai=fvvv_Evz_vz_vz();
                
                ai=ai+fvv_Evz_vzz()*omg1;
                
                ai=ai+A_times_kron_Q1_Qk(obj.debug,fplus_lmbd*obj.Tzz{rtp1},...
                    Jz_,Jzz_)*omg1;
                
                function fe=fvvv_Evz_vz_vz()
                    
                    fvvv_=obj.fvvv{rt,rtp1}(:,obj.expand.v{3});
                    
                    fe=A_times_kron_Q1_Qk(obj.debug,fvvv_,a0z_,a0z_,a0z_);
                                        
                    if any(obj.EU{3}(:))
                        
                        fe=fe+A_times_kron_Q1_Qk(obj.debug,fvvv_,a1z_,a1z_,...
                            a1z_)*obj.EU{3};
                                                
                    end
                    
                    fe=fe+fvvv_a0z_a1z_a1z_U2();
                    
                    function fe=fvvv_a0z_a1z_a1z_U2()
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            ];
                        
                        orders={[2,1,3],[2,3,1]};%[1,2,3], First added automatically
                        
                        options=[];
                        
                        fe=kron_Q1_Qk_times_A(obj.debug,a1z_,a1z_,obj.EU{2});
                        
                        fe=kron(a0z_,fe);
                        
                        fe=fvvv_*obj.sp(fe,matsizes,options,orders{:});
                        
                    end
                    
                end
                
                function fe=fvv_Evz_vzz()
                    
                    fvv_=obj.fvv{rt,rtp1}(:,obj.expand.v{2});
                    
                    EU2=obj.EU{2};
                    
                    fe=A_times_kron_Q1_Qk(obj.debug,fvv_,a0z_,a0zz_);
                    
                    if any(obj.EU{3}(:))
                        
                        fe=fe+A_times_kron_Q1_Qk(obj.debug,fvv_,a1z_,a1zz_)*obj.EU{3};
                        
                    end
                    
                    matsizes=[
                        obj.n.z,obj.n.z
                        obj.n.z,obj.n.z
                        obj.n.z,obj.n.z
                        ];
                    
                    orders={[1,3,2]};%[1,2,3], First added automatically
                    
                    options=[];
                    
                    P_kron_EU2_Jz=obj.sp(kron(EU2,Jz_),matsizes,options,orders{:});
                    % P(kron(EU2,Jz_))-kron(Jz_,EU2)
                    
                    fe=fe+A_times_kron_Q1_Qk(obj.debug,fvv_,a1z_,a1zz_)*P_kron_EU2_Jz;
                    
                    fe=fe+A_times_kron_Q1_Qk(obj.debug,fvv_,a0z_,a1zz_*EU2);
                    
                end
                
            end
            
            function ci=c_iterate()
                
                matsizes=[
                    obj.n.z,obj.n.z
                    obj.n.z,obj.n.z
                    obj.n.z,obj.n.z
                    ];
                
                orders={[2,1,3],[2,3,1]};%[1,2,3], First added automatically
                
                options=[];
                
                ci=obj.sp(kron(Jz_,obj.EU{2}),matsizes,options,orders{:})+...
                    obj.EU{3}+obj.kronall(Jz_,Jz_,Jz_);
                
                ci=W*ci*K;
                
            end
            
        end
        
        function obj=solve_fourth_order(obj)
            
            pb=obj.post.pb;
            
            bf=obj.post.bf;
            
            Jzzz_bulk=sparse(obj.n.e*(obj.k+1)+1,obj.n.z^3);
            
            a1zzz_bulk=sparse(obj.n.t+obj.n.pb+obj.n.e+obj.sig_in_v,obj.n.z^3);
            
            a0zzz_bulk=sparse(obj.n.pb+obj.n.e+obj.sig_in_v,obj.n.z^3);
            
            omg1=utils.cr.omega(obj.n.z,1);
            
            omg2=utils.cr.omega(obj.n.z,2);
            
            omg3=utils.cr.omega(obj.n.z,3);
            
            omg4=utils.cr.omega(obj.n.z,4);
            
            A=repmat({sparse(0)},1,obj.h);
            
            CBAR=cell(obj.h);
            
            K=obj.Comp_mat.z{4};
            
            W=obj.UnComp_mat.z{4};
            
            for rt=1:obj.h
                
                for rtp1=1:obj.h
                    
                    obj.Jzzz{rt,rtp1}=do_Jzzz();
                    
                    obj.a0zzz{rt,rtp1}=do_a0zzz();
                    
                    obj.a1zzz{rt,rtp1}=do_a1zzz();
                    
                    %--------------------------------
                    Jz_=obj.Jz{rt,rtp1};
                    
                    Jzz_=obj.Jzz{rt,rtp1};
                    
                    Jzzz_=obj.Jzzz{rt,rtp1};
                    % #
                    a0z_=obj.a0z{rt,rtp1};
                    
                    a0zz_=obj.a0zz{rt,rtp1};
                    
                    % a0zzz_=obj.a0zzz{rt,rtp1};
                    % #
                    a1z_=obj.a1z{rt,rtp1};
                    
                    a1zz_=obj.a1zz{rt,rtp1};
                    
                    a1zzz_=obj.a1zzz{rt,rtp1};
                    %--------------------------------
                    
                    A{rt}=A{rt}+a_iterate();
                    
                    CBAR{rt,rtp1}=c_iterate();
                    
                end
                % compress
                A{rt}=A{rt}*K;
                
            end
            
            obj.Tzzzz=solve_generalized_sylvester(A,obj.B,CBAR,...
                obj.fplus_lambda_bf,W,obj.gensylv_solve);
            
            function Jzzz=do_Jzzz()
                
                if rtp1==1
                    
                    tmp=obj.Tzzz{rt};
                    
                    Jzzz=[
                        tmp(pb,:)
                        Jzzz_bulk
                        ];
                    
                else
                    
                    Jzzz=obj.Jzzz{rt,1};
                    
                end
                
            end
            
            function a0zzz=do_a0zzz()
                
                Jz_=obj.Jz{rt,rtp1};
                
                Jzz_=obj.Jzz{rt,rtp1};
                
                Jzzz_=obj.Jzzz{rt,rtp1};
                
                Tz_=obj.Tz{rtp1};
                
                Tzz_=obj.Tzz{rtp1};
                
                Tzzz_=obj.Tzzz{rtp1};
                
                TJg=obj.cr3(Tzzz_,Tzz_,Tz_,Jzzz_,Jzz_,Jz_);
                
                a0zzz=[
                    TJg(bf,:)
                    obj.Tzzz{rt}
                    a0zzz_bulk
                    ];
            end
            
            function a1zzz=do_a1zzz()
                
                a1zzz=[obj.Tzzz{rtp1}(bf,:)
                    a1zzz_bulk];
                
            end
            
            function ai=a_iterate()
                
                fplus_lmbd=obj.fv{rt,rtp1}(:,obj.posv.bf)*obj.lambda.bf;
                
                fvvv_=obj.fvvv{rt,rtp1}(:,obj.expand.v{3});
                
                fvv_=obj.fvv{rt,rtp1}(:,obj.expand.v{2});
                
                ai=fvvvv_Evz_vz_vz_vz();
                
                ai=ai+fvvv_*Evz2_vzz()*omg2;
                
                ai=ai+fvv_*Evz_vzzz()*omg3;
                
                ai=ai+fvv_*Evzz_vzz()*omg4;
                
                ai=ai+fplus_lmbd*obj.Tzzz{rtp1}*(kron(obj.EU{2},Jzz_)*omg2+...
                    obj.kronall(Jz_,Jz_,Jzz_)*omg2);
                
                ai=ai+fplus_lmbd*obj.Tzz{rtp1}*(kron(Jzz_,Jzz_)*omg4+...
                    kron(Jz_,Jzzz_)*omg3);
                
                function fe=fvvvv_Evz_vz_vz_vz()
                    
                    fvvvv_=obj.fvvvv{rt,rtp1}(:,obj.expand.v{4});
                    
                    EU2=obj.EU{2};
                    
                    a1z3=obj.kronall(a1z_,a1z_,a1z_);
                    
                    fe=A_times_kron_Q1_Qk(obj.debug,fvvvv_,a0z_,a0z_,a0z_,a0z_);
                    
                    fe=fe+A_times_kron_Q1_Qk(obj.debug,fvvvv_,a1z_,a1z3)*obj.EU{4};
                    
                    options=[];
                    
                    matsizes=[
                        obj.n.v,obj.n.z
                        obj.n.v,obj.n.z
                        obj.n.v,obj.n.z
                        obj.n.v,obj.n.z
                        ];
                    
                    e=do_a0z_a0z_a1z_a1z_EU2();
                    
                    e=e+do_a0z_a1zU_a1zU_a1zU();
                    
                    fe=fe+fvvvv_*e;
                    
                    %----------------------------------
                    function p=do_a0z_a1zU_a1zU_a1zU()
                        
                        p=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        orders={[2,1,3,4],[2,3,1,4],[2,3,4,1]};%[1,2,3,4], First added automatically
                        
                        p=obj.sp(kron(a0z_,a1z3*obj.EU{3}),matsizes,options,orders{:});
                        
                    end
                    
                    function p=do_a0z_a0z_a1z_a1z_EU2()
                        
                        orders={[1,3,2,4],[1,3,4,2],[3,1,4,2],[3,4,1,2]};%[1,2,3,4], First added automatically
                        
                        p=obj.sp(obj.kronall(a0z_,a0z_,kron(a1z_,a1z_)*EU2),...
                            matsizes,options,orders{:});
                        
                    end
                    %----------------------------------
                    
                end
                
                function e=Evz2_vzz()
                    
                    e=obj.kronall(a0z_,a0z_,a0zz_);
                    
                    e=e+kron(kron(a1z_,a1z_)*obj.EU{2},a0zz_);
                    
                    e=e+kron(speye(obj.n.v^2),a1zz_)*(do_a0z_a1zU_U_Jz()+...
                        do_a1z2_U2_U_Jz());
                    
                    e=e+obj.kronall(a0z_,a0z_,a1zz_*obj.EU{2});
                    
                    e=e+do_a0z_a1zU_a1zz_U2();
                    
                    e=e+obj.kronall(a1z_,a1z_,a1zz_)*obj.EU{4};
                    
                    function p=do_a1z2_U2_U_Jz()
                        
                        p=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        p=kron(obj.kronall(a1z_,a1z_,speye(obj.n.z))*obj.EU{3},Jz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.v^2,obj.n.z^2
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,3,2]};%[1,2,3], First added automatically
                        
                        p=obj.sp(p,matsizes,options,orders{:});
                        
                    end
                    
                    function p=do_a0z_a1zU_U_Jz()
                        
                        p=obj.kronall(a0z_,kron(a1z_,speye(obj.n.z))*obj.EU{2},Jz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[2,1,3,4],[2,1,4,3],[1,2,4,3]};%[1,2,3,4], First added automatically
                        
                        p=obj.sp(p,matsizes,options,orders{:});
                        
                    end
                    
                    function p=do_a0z_a1zU_a1zz_U2()
                        
                        p=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z^2
                            ];
                        
                        orders={[2,1,3]};%[1,2,3], First added automatically
                        
                        p=obj.sp(kron(a0z_,kron(a1z_,a1zz_)*obj.EU{3}),matsizes,options,orders{:});
                        
                    end
                    
                end
                
                function e=Evz_vzzz()
                    
                    e=kron(a0z_,a1zzz_*do_U2_Jz());
                    
                    if any(obj.EU{3}(:))
                        
                        e=e+kron(a0z_,a1zzz_*obj.EU{3});
                        
                    end
                    
                    e=e+do_a1zU_a1zzz_U_Jz2();
                    
                    e=e+do_a1zU_a1zzz_U2_Jz();
                    
                    e=e+kron(a1z_,a1zzz_)*obj.EU{4};
                    
                    e=e+kron(a1z_,a1zz_)*kron(obj.EU{2},Jzz_)*kron(speye(obj.n.z),omg1);
                    
                    function p=do_a1zU_a1zzz_U2_Jz()
                        
                        p=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        p=obj.kronall(obj.EU{3},Jz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,2,4,3],[1,4,2,3]};%[1,2,3,4], First added automatically
                        
                        p=kron(a1z_,a1zzz_)*obj.sp(p,matsizes,options,orders{:});
                        
                    end
                    
                    function p=do_a1zU_a1zzz_U_Jz2()
                        
                        p=obj.kronall(obj.EU{2},Jz_,Jz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,3,2,4],[1,3,4,2]};%[1,2,3,4], First added automatically
                        
                        p=kron(a1z_,a1zzz_)*obj.sp(p,matsizes,options,orders{:});
                        
                    end
                    
                    function p=do_U2_Jz()
                        
                        p=kron(obj.EU{2},Jz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,3,2],[3,1,2]};%[1,2,3], First added automatically
                        
                        p=obj.sp(p,matsizes,options,orders{:});
                        
                    end
                    
                end
                
                function e=Evzz_vzz()
                    
                    e=kron(a0zz_,a0zz_);
                    
                    e=e+do_a0zz_a1zzU2();
                    
                    e=e+do_a1zz_Jz_U_power2();
                    
                    e=e+do_a1zz_Jz_U_a1zz_U2();
                    
                    e=e+kron(a1zz_,a1zz_)*obj.EU{4};
                    
                    function p=do_a1zz_Jz_U_a1zz_U2()
                        
                        p=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        p=kron(Jz_,obj.EU{3});
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z^2,obj.n.z^2
                            ];
                        
                        orders={[2,1,3]};%[1,2,3], First added automatically
                        
                        % inner part
                        %-------------
                        p=kron(a1zz_,a1zz_)*obj.sp(p,matsizes,options,orders{:});
                        
                        % outer part
                        %-------------
                        matsizes=[
                            obj.n.v,obj.n.z^2
                            obj.n.v,obj.n.z^2
                            ];
                        
                        orders={[2,1]};%[1,2], First added automatically
                        
                        p=obj.sp(p,matsizes,options,orders{:});
                        
                    end
                    
                    function p=do_a1zz_Jz_U_power2()
                        
                        p=obj.kronall(Jz_,obj.EU{2},Jz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[2,1,3,4],[2,1,4,3],[1,2,4,3]};%[1,2,3,4], First added automatically
                        
                        p=kron(a1zz_,a1zz_)*obj.sp(p,matsizes,options,orders{:});
                        
                    end
                    
                    function p=do_a0zz_a1zzU2()
                        
                        p=kron(a0zz_,a1zz_*obj.EU{2});
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z^2
                            obj.n.v,obj.n.z^2
                            ];
                        
                        orders={[2,1]};%[1,2], First added automatically
                        
                        p=obj.sp(p,matsizes,options,orders{:});
                        
                    end
                    
                end
                
            end
            
            function ci=c_iterate()
                
                options=[];
                
                matsizes=[
                    obj.n.z,obj.n.z
                    obj.n.z,obj.n.z
                    obj.n.z,obj.n.z
                    obj.n.z,obj.n.z
                    ];
                
                
                ci=obj.kronall(Jz_,Jz_,Jz_,Jz_)+do_Jz_Jz_EU2()+do_Jz_EU3()+obj.EU{4};
                
                ci=W*ci*K;
                
                function p=do_Jz_Jz_EU2()
                    
                    orders={[1,3,2,4],[1,3,4,2],[3,1,4,2],[3,4,1,2]};%[1,2,3,4], First added automatically
                    
                    p=obj.sp(obj.kronall(Jz_,Jz_,obj.EU{2}),matsizes,options,orders{:});
                    
                end
                
                function p=do_Jz_EU3()
                    
                    p=sparse(0);
                    
                    if ~any(obj.EU{3}(:))
                        
                        return
                        
                    end
                    
                    orders={[2,1,3,4],[2,3,1,4],[2,3,4,1]};%[1,2,3,4], First added automatically
                    
                    p=obj.sp(kron(Jz_,obj.EU{3}),matsizes,options,orders{:});
                    
                end
                
            end
            
        end
        
        function obj=solve_fifth_order(obj)
            
            pb=obj.post.pb;
            
            bf=obj.post.bf;
            
            Jzzzz_bulk=sparse(obj.n.e*(obj.k+1)+1,obj.n.z^4);
            
            a1zzzz_bulk=sparse(obj.n.t+obj.n.pb+obj.n.e+obj.sig_in_v,obj.n.z^4);
            
            a0zzzz_bulk=sparse(obj.n.pb+obj.n.e+obj.sig_in_v,obj.n.z^4);
            
            omg1=utils.cr.omega(obj.n.z,1);
            
            omg2=utils.cr.omega(obj.n.z,2);
            
            omg3=utils.cr.omega(obj.n.z,3);
            
            omg4=utils.cr.omega(obj.n.z,4);
            
            omg5=utils.cr.omega(obj.n.z,5);
            
            omg6=utils.cr.omega(obj.n.z,6);
            
            omg7=utils.cr.omega(obj.n.z,7);
            
            omg8=utils.cr.omega(obj.n.z,8);
            
            omg9=utils.cr.omega(obj.n.z,9);
            
            A=repmat({sparse(0)},1,obj.h);
            
            CBAR=cell(obj.h);
            
            K=obj.Comp_mat.z{5};
            
            W=obj.UnComp_mat.z{5};
            
            for rt=1:obj.h
                
                for rtp1=1:obj.h
                    
                    obj.Jzzzz{rt,rtp1}=do_Jzzzz();
                    
                    Jz_=obj.Jz{rt,rtp1};
                    
                    Jzz_=obj.Jzz{rt,rtp1};
                    
                    Jzzz_=obj.Jzzz{rt,rtp1};
                    
                    Jzzzz_=obj.Jzzzz{rt,rtp1};
                    
                    obj.a0zzzz{rt,rtp1}=do_a0zzzz();
                    
                    obj.a1zzzz{rt,rtp1}=do_a1zzzz();
                    
                    %--------------------------------
                    % #
                    a0z_=obj.a0z{rt,rtp1};
                    
                    a0zz_=obj.a0zz{rt,rtp1};
                    
                    a0zzz_=obj.a0zzz{rt,rtp1};
                    
                    a0zzzz_=obj.a0zzzz{rt,rtp1};
                    % #
                    a1z_=obj.a1z{rt,rtp1};
                    
                    a1zz_=obj.a1zz{rt,rtp1};
                    
                    a1zzz_=obj.a1zzz{rt,rtp1};
                    
                    a1zzzz_=obj.a1zzzz{rt,rtp1};
                    %--------------------------------
                    
                    A{rt}=A{rt}+a_iterate();
                    
                    CBAR{rt,rtp1}=c_iterate();
                    
                end
                % compress
                A{rt}=A{rt}*K;
                
            end
            
            obj.Tzzzzz=solve_generalized_sylvester(A,obj.B,CBAR,...
                obj.fplus_lambda_bf,W,obj.gensylv_solve);
            
            function Jzzzz=do_Jzzzz()
                
                if rtp1==1
                    
                    Jzzzz=[
                        obj.Tzzzz{rt}(pb,:)
                        Jzzzz_bulk
                        ];
                    
                else
                    
                    Jzzzz=obj.Jzzzz{rt,1};
                    
                end
                
            end
            
            function a0zzzz=do_a0zzzz()
                
                Tz_=obj.Tz{rtp1};
                
                Tzz_=obj.Tzz{rtp1};
                
                Tzzz_=obj.Tzzz{rtp1};
                
                Tzzzz_=obj.Tzzzz{rtp1};
                
                if obj.use_kron
                    
                    TJg=Tzzzz_*obj.kronall(Jz_,Jz_,Jz_,Jz_)+...
                        Tzzz_*obj.kronall(Jz_,Jz_,Jzz_)*omg2+...
                        Tzz_*kron(Jz_,Jzzz_)*omg3+...
                        Tzz_*kron(Jzz_,Jzz_)*omg4+...
                        Tz_*Jzzzz_;
                    
                else
                    
                    TJg=obj.cr4(Tzzzz_,Tzzz_,Tzz_,Tz_,Jzzzz_,Jzzz_,Jzz_,Jz_);
                    
                end
                
                a0zzzz=[
                    TJg(bf,:)
                    obj.Tzzzz{rt}
                    a0zzzz_bulk
                    ];
            end
            
            function a1zzzz=do_a1zzzz()
                
                a1zzzz=[obj.Tzzzz{rtp1}(bf,:)
                    a1zzzz_bulk];
                
            end
            
            function ai=a_iterate()
                
                Tz_=obj.Tz{rtp1};
                
                Tzz_=obj.Tzz{rtp1};
                
                Tzzz_=obj.Tzzz{rtp1};
                
                Tzzzz_=obj.Tzzzz{rtp1};
                
                ai=obj.fv{rt,rtp1}(:,obj.posv.bf)*obj.lambda.bf*...
                    obj.cr5([],Tzzzz_,Tzzz_,Tzz_,Tz_,...
                    [],Jzzzz_,Jzzz_,Jzz_,Jz_);
                
                ai=ai+...
                    fvvvvv_Evz_vz_vz_vz_vz()+...
                    fvvvv_Evz3_vzz()+...
                    fvvv_Evz_vz_vzzz()+...
                    fvvv_Evz_vzz_vzz()+...
                    fvv_Evz_vzzzz()+...
                    fvv_Evzz_vzzz()+...
                    fv_a1zzzz_Jz_U2_U3_Jzz_a1zzz_U2_Jzzz();
                
                function fe=fvvvvv_Evz_vz_vz_vz_vz()
                    
                    fvvvvv_=obj.fvvvvv{rt,rtp1}(:,obj.expand.v{5});
                    
                    fe=A_times_kron_Q1_Qk(obj.debug,fvvvvv_,a0z_,a0z_,a0z_,a0z_,a0z_);
                    
                    fe=fe+fvvvvv_a0z_a1z4_U4();
                                        
                    e=a0z2_a1z3U3();
                    
                    etmp=a0z3_a1z2U2();
                    
                    e=e+etmp;
                    
                    if any(obj.EU{5}(:))
                        
                        e=e+obj.kronall(a1z_,a1z_,a1z_,a1z_,a1z_)*obj.EU{5};
                        
                    end
                    
                    fe=fe+fvvvvv_*e;
                    
                    function fe=fvvvvv_a0z_a1z4_U4()
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v^4,obj.n.z^4
                            ];
                        
                        orders={[2,1]};%[1,2], First added automatically
                        
                        % e=obj.kronall(a1z_,a1z_,a1z_,a1z_)*obj.EU{4};
                        e=kron_Q1_Qk_times_A(obj.debug,a1z_,a1z_,a1z_,a1z_,obj.EU{4});
                        
                        e=kron(a0z_,e);
                        
                        e=obj.sp(e,matsizes,options,orders{:});
                        
                        fe=fvvvvv_*e;
                        
                    end
                    
                    function e=a0z2_a1z3U3()
                        
                        e=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        e=obj.kronall(a0z_,a0z_,obj.kronall(a1z_,a1z_,a1z_)*obj.EU{3});
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.v^3,obj.n.z^3
                            ];
                        
                        orders={[1,3,2],[3,1,2]};%[1,2,3], First added automatically
                        
                        e=obj.sp(e,matsizes,options,orders{:});
                        
                    end
                    
                    function e=a0z3_a1z2U2()
                        
                        e=obj.kronall(a0z_,a0z_,a0z_,...
                            kron_Q1_Qk_times_A(obj.debug,a1z_,a1z_,obj.EU{2}));
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.v^2,obj.n.z^2
                            ];
                        
                        orders={[1,2,4,3],[1,4,2,3],[4,1,2,3]};%[1,2,3,4], First added automatically
                        
                        e=obj.sp(e,matsizes,options,orders{:});
                        
                    end
                    
                end
                
                function fe=fvvvv_Evz3_vzz()
                    
                    fvvvv_=obj.fvvvv{rt,rtp1}(:,obj.expand.v{4});
                                        
                    fe=A_times_kron_Q1_Qk(obj.debug,fvvvv_,a0z_,a0z_,a0z_,a0zz_+a1zz_*obj.EU{2});
                    
                    e=I_a1zz_a0z2_a1zU_U_Jz_or_U2();
                    
                    e=e+a0z_a1zU_a1zU_a0zz();
                    
                    e=e+a0z_a1zU_a1zU_a1zz_U_Jz();
                    
                    e=e+a0z_a1zU_a1zU_a1zzU2();
                    
                    e=e+a1zU_a1zU_a1zU_a0zz();
                    
                    e=e+a1zU_a1zU_a1zU_a1zz_U_Jz();
                    
                    e=e+a1zU_a1zU_a1zU_a1zz_U2();
                    
                    fe=fe+fvvvv_*e*omg5;
                    
                    function p=I_a1zz_a0z2_a1zU_U_Jz_or_U2()
                        
                        p1=do_U2_Jz();
                        
                        p2=do_U3();
                        
                        p=kron(speye(obj.n.v^3),a1zz_)*(p1+p2);
                        
                        function p=do_U3()
                            
                            p=sparse(0);
                            
                            if ~any(obj.EU{3}(:))
                                
                                return
                                
                            end
                            
                            p=obj.kronall(a0z_,a0z_,kron(a1z_,speye(obj.n.z^2))*obj.EU{3});
                            
                            options=[];
                            
                            matsizes=[
                                obj.n.v,obj.n.z
                                obj.n.v,obj.n.z
                                obj.n.v,obj.n.z
                                obj.n.z^2,obj.n.z^2
                                ];
                            
                            orders={[1,3,2,4],[3,1,2,4]};%[1,2,3,4], First added automatically
                            
                            p=obj.sp(p,matsizes,options,orders{:});
                            
                        end
                        
                        function p=do_U2_Jz()
                            
                            p=obj.kronall(a0z_,a0z_,kron(a1z_,speye(obj.n.z^2))*kron(obj.EU{2},Jz_));
                            
                            options=[];
                            
                            matsizes=[
                                obj.n.v,obj.n.z
                                obj.n.v,obj.n.z
                                obj.n.v,obj.n.z
                                obj.n.z,obj.n.z
                                obj.n.z,obj.n.z
                                ];
                            
                            orders={[1,3,2,4,5],[3,1,2,4,5],[1,2,3,5,4],...
                                [1,3,2,5,4],[3,1,2,5,4]};%[1,2,3,4,5], First added automatically
                            
                            p=obj.sp(p,matsizes,options,orders{:});
                            
                        end
                        
                    end
                    
                    function p=a0z_a1zU_a1zU_a0zz()
                        
                        p=kron(a0z_,kron(a1z_,a1z_)*obj.EU{2});
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            ];
                        
                        orders={[2,1,3],[2,3,1]};%[1,2,3], First added automatically
                        
                        p=obj.sp(p,matsizes,options,orders{:});
                        
                        p=kron(p,a0zz_);
                        
                    end
                    
                    function p=a0z_a1zU_a1zU_a1zz_U_Jz()
                        
                        p=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        p=obj.kronall(a1z_,a1z_,eye(obj.n.z))*obj.EU{3};
                        
                        p=obj.kronall(a0z_,p,Jz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[2,1,3,4,5],[2,3,1,4,5],[2,1,3,5,4],...
                            [2,3,1,5,4],[1,2,3,5,4]};%[1,2,3,4,5], First added automatically
                        
                        p=obj.sp(p,matsizes,options,orders{:});
                        
                        p=kron(eye(obj.n.v^3),a1zz_)*p;
                        
                    end
                    
                    function p=a0z_a1zU_a1zU_a1zzU2()
                        
                        p=obj.kronall(a1z_,a1z_,a1zz_)*obj.EU{4};
                        
                        p=kron(a0z_,p);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z^2
                            ];
                        
                        orders={[2,1,3,4],[2,3,1,4]};%[1,2,3,4], First added automatically
                        
                        p=obj.sp(p,matsizes,options,orders{:});
                        
                    end
                    
                    function p=a1zU_a1zU_a1zU_a0zz()
                        
                        p=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        p=obj.kronall(a1z_,a1z_,a1z_)*obj.EU{3};
                        
                        p=kron(p,a0zz_);
                        
                    end
                    
                    function p=a1zU_a1zU_a1zU_a1zz_U_Jz()
                        
                        p=kron(obj.EU{4},Jz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.z^3,obj.n.z^3
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,3,2]};%[1,2,3], First added automatically
                        
                        p=obj.sp(p,matsizes,options,orders{:});
                        
                        p=obj.kronall(a1z_,a1z_,a1z_,a1zz_)*p;
                        
                    end
                    
                    function p=a1zU_a1zU_a1zU_a1zz_U2()
                        
                        p=sparse(0);
                        
                        if ~any(obj.EU{5}(:))
                            
                            return
                            
                        end
                        
                        p=obj.kronall(a1z_,a1z_,a1z_,a1zz_)*obj.EU{5};
                        
                    end
                    
                end
                
                function fe=fvvv_Evz_vz_vzzz()
                    
                    fvvv_=obj.fvvv{rt,rtp1}(:,obj.expand.v{3});
                    
                    fe=A_times_kron_Q1_Qk(obj.debug,fvvv_,a0z_,a0z_,a0zzz_);
                    
                    fe=fe+A_times_kron_Q1_Qk(obj.debug,fvvv_,a0z_,a0z_,do_a1zzz_U2_Jz());
                    
                    fe=fe+fvvv_a0z_a0z_a1zzzU3();
                    
                    e=a0z_a1zU_a1zzz_U_Jz2();
                    
                    e=e+a0z_a1zU_a1zzz_U2_Jz();
                    
                    e=e+a0z_a1zU_a1zzzU3();
                    
                    e=e+a0z_a1zU_a1zz_U_Jzz_omg1();
                    
                    fe=fe+A_times_kron_Q1_Qk(obj.debug,fvvv_,...
                        kron_Q1_Qk_times_A(obj.debug,a1z_,a1z_,obj.EU{2}),...
                        a0zzz_);
                    
                    e=e+a1zU_a1zU_a1zzz_U_Jz2();
                    
                    e=e+a1zU_a1zU_a1zzz_U2_Jz();
                    
                    if any(obj.EU{5}(:))
                        
                        fe=fe+A_times_kron_Q1_Qk(obj.debug,fvvv_,a1z_,a1z_,a1zzz_)*obj.EU{5};
                        
                    end
                    
                    fe=fe+fvvv_a1zU_a1zU_a1zz_U_Jzz_omg1();
                    
                    fe=fe+fvvv_*e;
                    
                    fe=fe*omg6;
                    
                    function f=fvvv_a1zU_a1zU_a1zz_U_Jzz_omg1()
                        
                        f=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        f=A_times_kron_Q1_Qk(obj.debug,fvvv_,a1z_,a1z_,a1zz_);
                        
                        f=A_times_kron_Q1_Qk(obj.debug,f,obj.EU{3},Jzz_);
                        
                        f=A_times_kron_Q1_Qk(obj.debug,f,speye(obj.n.z^2),omg1);
                        
                    end
                    
                    function e=a1zU_a1zU_a1zzz_U2_Jz()
                        
                        e=obj.kronall(a1z_,a1z_,speye(obj.n.z^2))*obj.EU{4};
                        
                        e=kron(e,Jz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.v^2,obj.n.z^2
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,2,4,3],[1,4,2,3]};%[1,2,3,4], First added automatically
                        
                        e=obj.sp(e,matsizes,options,orders{:});
                        
                        e=kron(speye(obj.n.v^2),a1zzz_)*e;
                        
                    end
                    
                    function e=a1zU_a1zU_a1zzz_U_Jz2()
                        
                        e=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        e=obj.kronall(a1z_,a1z_,speye(obj.n.z))*obj.EU{3};
                        
                        e=obj.kronall(e,Jz_,Jz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.v^2,obj.n.z^2
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,3,2,4],[1,3,4,2]};%[1,2,3,4], First added automatically
                        
                        e=obj.sp(e,matsizes,options,orders{:});
                        
                        e=kron(speye(obj.n.v^2),a1zzz_)*e;
                        
                    end
                    
                    function e=do_a1zzz_U2_Jz()
                        
                        e=kron(obj.EU{2},Jz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,3,2],[3,1,2]};%[1,2,3], First added automatically
                        
                        e=obj.sp(e,matsizes,options,orders{:});
                        
                        e=a1zzz_*e;
                        
                    end
                    
                    function f=fvvv_a0z_a0z_a1zzzU3()
                        
                        f=sparse(0);
                        
                        if any(obj.EU{3}(:))
                            
                            f=A_times_kron_Q1_Qk(obj.debug,fvvv_,a0z_,a0z_,a1zzz_*obj.EU{3});
                            
                        end
                        
                    end
                    
                    function e=a0z_a1zU_a1zzz_U_Jz2()
                        
                        e=kron_Q1_Qk_times_A(obj.debug,a1z_,speye(obj.n.z),obj.EU{2});
                        
                        e=obj.kronall(a0z_,e,Jz_,Jz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,2,4,3,5],[1,2,4,5,3],...
                            [2,1,3,4,5],[2,1,4,3,5],[2,1,4,5,3]};%[1,2,3,4,5], First added automatically
                        
                        e=obj.sp(e,matsizes,options,orders{:});
                        
                        e=kron_Q1_Qk_times_A(obj.debug,speye(obj.n.v^2),a1zzz_,e);
                        
                    end
                    
                    function e=a0z_a1zU_a1zzz_U2_Jz()
                        
                        e=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        e=kron_Q1_Qk_times_A(obj.debug,a1z_,speye(obj.n.z^2),obj.EU{3});
                        
                        e=obj.kronall(a0z_,e,Jz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,2,3,5,4],[1,2,5,3,4],...
                            [2,1,3,4,5],[2,1,3,5,4],[2,1,5,3,4]};%[1,2,3,4,5], First added automatically
                        
                        e=obj.sp(e,matsizes,options,orders{:});
                        
                        e=kron_Q1_Qk_times_A(obj.debug,speye(obj.n.v^2),a1zzz_,e);
                        
                    end
                    
                    function e=a0z_a1zU_a1zzzU3()
                        
                        e=kron(a0z_,kron_Q1_Qk_times_A(obj.debug,a1z_,a1zzz_,obj.EU{4}));
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z^3
                            ];
                        
                        orders={[2,1,3]};%[1,2,3], First added automatically
                        
                        e=obj.sp(e,matsizes,options,orders{:});
                        
                    end
                    
                    function e=a0z_a1zU_a1zz_U_Jzz_omg1()
                        
                        e=kron_Q1_Qk_times_A(obj.debug,a1z_,speye(obj.n.z),obj.EU{2});
                        
                        e=kron_Q1_Qk_times_A(obj.debug,a0z_,e,Jzz_,...
                            kron(speye(obj.n.z^2),omg1)...
                            );
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.z^2,obj.n.z^3
                            ];
                        
                        orders={[2,1,3]};%[1,2,3], First added automatically
                        
                        e=obj.sp(e,matsizes,options,orders{:});
                        
                        e=kron_Q1_Qk_times_A(obj.debug,speye(obj.n.v^2),a1zz_,e);
                        
                    end
                    
                end
                
                function fe=fvvv_Evz_vzz_vzz()
                    
                    fvvv_=obj.fvvv{rt,rtp1}(:,obj.expand.v{3});
                    
                    fe=A_times_kron_Q1_Qk(obj.debug,fvvv_,a0z_,a0zz_,a0zz_);
                    
                    e=a0z_a0zz_a1zzU2();
                    
                    e=e+a0z_a1zz_Jz_U_a1zzU2();
                    
                    e=e+a0z_a1zz_Jz_U_a1zz_U_Jz();
                    
                    fe=fe+A_times_kron_Q1_Qk(obj.debug,fvvv_,a0z_,...
                        kron_Q1_Qk_times_A(obj.debug,a1zz_,a1zz_,obj.EU{4})...
                        );
                    
                    e=e+a1zU_a1zz_U_Jz_a0zz();
                    
                    e=e+a1zU_a1zzU2_a0zz();
                    
                    e=e+a1zU_a1zzU2_a1zz_U_Jz();
                    
                    e=e+a1zU_a1zz_U_Jz_a1zz_U_Jz();
                    
                    if any(obj.EU{5}(:))
                        
                        fe=fe+A_times_kron_Q1_Qk(obj.debug,fvvv_,a1z_,a1zz_,...
                            a1zz_)*obj.EU{5};
                        
                    end
                    
                    fe=fe+fvvv_*e;
                    
                    fe=fe*omg7;
                    
                    function e=a1zU_a1zz_U_Jz_a1zz_U_Jz()
                        
                        e=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        e=obj.kronall(obj.EU{3},Jz_,Jz_);
                        % put in original form
                        %---------------------
                        options=struct('skip_first',true);
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,4,2,5,3]};%[1,2,3,4,5], First NOT added
                        
                        e=obj.sp(e,matsizes,options,orders{:});
                        
                        % Now order according to permutations
                        %-------------------------------------
                        orders={[1,2,3,4,5],[1,2,3,5,4],...
                            [1,3,2,4,5],[1,3,2,5,4]};%[1,2,3,4,5], First NOT added automatically
                        
                        e=obj.sp(e,matsizes,options,orders{:});
                        
                        % then multiply with coefficient
                        %-------------------------------
                        e=obj.kronall(a1z_,a1zz_,a1zz_)*e;
                        
                    end
                    
                    function e=a1zU_a1zzU2_a1zz_U_Jz()
                        
                        e=kron(obj.EU{4},Jz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.z^3,obj.n.z^3
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,3,2]};%[1,2,3], First added automatically
                        
                        e=obj.sp(e,matsizes,options,orders{:});
                        
                        e=obj.kronall(a1z_,a1zz_,a1zz_)*e;
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z^2
                            obj.n.v,obj.n.z^2
                            ];
                        
                        orders={[1,3,2]};%[1,2,3], First added automatically
                        
                        e=obj.sp(e,matsizes,options,orders{:});
                        
                        
                    end
                    
                    function e=a1zU_a1zzU2_a0zz()
                        
                        e=kron(a1z_,a1zz_)*obj.EU{3};
                        
                        e=kron(e,a0zz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z^2
                            obj.n.v,obj.n.z^2
                            ];
                        
                        orders={[1,3,2]};%[1,2,3], First added automatically
                        
                        e=obj.sp(e,matsizes,options,orders{:});
                        
                    end
                    
                    function e=a1zU_a1zz_U_Jz_a0zz()
                        
                        e=kron(obj.EU{2},Jz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,3,2]};%[1,2,3], First added automatically
                        
                        e=obj.sp(e,matsizes,options,orders{:});
                        
                        e=kron(a1z_,a1zz_)*e;
                        
                        e=kron(e,a0zz_);
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z^2
                            obj.n.v,obj.n.z^2
                            ];
                        
                        orders={[1,3,2]};%[1,2,3], First added automatically
                        
                        e=obj.sp(e,matsizes,options,orders{:});
                        
                    end
                    
                    function e=a0z_a1zz_Jz_U_a1zz_U_Jz()
                        
                        e=obj.kronall(Jz_,obj.EU{2},Jz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,2,4,3],[2,1,3,4],[2,1,4,3]};%[1,2,3,4], First added automatically
                        
                        e=obj.sp(e,matsizes,options,orders{:});
                        
                        e=kron(a1zz_,a1zz_)*e;
                        
                        e=kron(a0z_,e);
                        
                        
                    end
                    
                    function e=a0z_a1zz_Jz_U_a1zzU2()
                        
                        e=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        e=kron(Jz_,obj.EU{3});
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z^2,obj.n.z^2
                            ];
                        
                        orders={[2,1,3]};%[1,2,3], First added automatically
                        
                        e=obj.sp(e,matsizes,options,orders{:});
                        
                        e=kron(a1zz_,a1zz_)*e;
                        
                        matsizes=[
                            obj.n.v,obj.n.z^2
                            obj.n.v,obj.n.z^2
                            ];
                        
                        orders={[2,1]};%[1,2], First added automatically
                        
                        e=obj.sp(e,matsizes,options,orders{:});
                        
                        e=kron(a0z_,e);
                        
                    end
                    
                    function e=a0z_a0zz_a1zzU2()
                        
                        e=obj.kronall(a0z_,a0zz_,a1zz_*obj.EU{2});
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z^2
                            obj.n.v,obj.n.z^2
                            ];
                        
                        orders={[1,3,2]};%[1,2,3], First added automatically
                        
                        e=obj.sp(e,matsizes,options,orders{:});
                        
                    end
                    
                end
                
                function fe=fvv_Evz_vzzzz()
                    
                    e=kron(a0z_,a0zzzz_);
                    
                    e=e+kron(a0z_,a1zzzz_*...
                        (a1zzzz_Jz2_U2()+Jz_U3()+obj.EU{4})+...
                        a1zzz_*kron(obj.EU{2},Jzz_)*omg2);
                    
                    e=e+kron(a1z_,a1zzzz_)*(...
                        U_U_Jz3()+U_U2_Jz2()+U_U3_Jz()+obj.EU{5}...
                        );
                    
                    e=e+kron(a1z_,a1zzz_)*U_U_Jz_U2_Jzz();
                    
                    e=e+kron(a1z_,a1zz_)*kron(obj.EU{2},Jzzz_)*kron(speye(obj.n.z),omg3);
                    
                    fe=obj.fvv{rt,rtp1}(:,obj.expand.v{2})*e*omg8;
                    
                    function p=U_U_Jz_U2_Jzz()
                        % This needs to be redone.
                        
                        p=kron(obj.EU{2},Jz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,3,2],};%[1,2,3], First added automatically
                        
                        p=obj.sp(p,matsizes,options,orders{:});
                        
                        p=p+obj.EU{3};
                        
                        p=kron(p,Jzz_)*kron(speye(obj.n.z),omg2);
                        
                    end
                    
                    function p=U_U3_Jz()
                        
                        p=kron(obj.EU{4},Jz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,2,3,5,4],[1,2,5,3,4],[1,5,2,3,4]};%[1,2,3,4,5], First added automatically
                        
                        p=obj.sp(p,matsizes,options,orders{:});
                        
                    end
                    
                    function p=U_U2_Jz2()
                        
                        p=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        p=obj.kronall(obj.EU{3},Jz_,Jz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,2,4,3,5],[1,4,2,3,5],[4,1,2,3,5],...
                            [1,2,4,5,3],[1,4,5,2,3]};%[1,2,3,4,5], First added automatically
                        
                        p=obj.sp(p,matsizes,options,orders{:});
                        
                        
                    end
                    
                    function p=U_U_Jz3()
                        
                        p=obj.kronall(obj.EU{2},Jz_,Jz_,Jz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,3,2,4,5],[1,3,4,2,5],[1,3,4,5,2]};%[1,2,3,4,5], First added automatically
                        
                        p=obj.sp(p,matsizes,options,orders{:});
                        
                    end
                    
                    function p=Jz_U3()
                        
                        p=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        p=kron(Jz_,obj.EU{3});
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[2,1,3,4],[2,3,1,4],[2,3,4,1]};%[1,2,3,4], First added automatically
                        
                        p=obj.sp(p,matsizes,options,orders{:});
                        
                    end
                    
                    function p=a1zzzz_Jz2_U2()
                        
                        p=obj.kronall(Jz_,Jz_,obj.EU{2});
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,3,2,4],[1,3,4,2],[3,1,4,2],[3,4,1,2]};%[1,2,3,4], First added automatically
                        
                        p=obj.sp(p,matsizes,options,orders{:});
                        
                    end
                    
                end
                
                function fe=fvv_Evzz_vzzz()
                    
                    e=kron(a0zz_,a0zzz_);
                    
                    e=e+kron(a0zz_,a1zzz_*(U2_Jz()+obj.EU{3}));
                    
                    e=e+kron(a1zz_,a1zzz_)*(Jz_U_U_Jz2()+Jz_U_U2_Jz()+Jz_U_U3());
                    
                    e=e+kron(a1zz_,a1zz_)*Jz_U_U_Jzz()*kron(speye(obj.n.z^2),omg1);
                    
                    e=e+kron(a1zz_*obj.EU{2},a0zzz_);
                    
                    e=e+kron(a1zz_,a1zzz_)*(U2_U_Jz2()+U2_U2_Jz()+obj.EU{5});
                    
                    if any(obj.EU{3}(:))
                        
                        e=e+kron(a1zz_,a1zz_)*kron(obj.EU{3},Jzz_)*kron(speye(obj.n.z^2),omg1);
                        
                    end
                    
                    fe=obj.fvv{rt,rtp1}(:,obj.expand.v{2})*e*omg9;
                    
                    function p=U2_U2_Jz()
                        
                        p=kron(obj.EU{4},Jz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.z^2,obj.n.z^2
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,2,4,3],[1,4,2,3]};%[1,2,3,4], First added automatically
                        
                        p=obj.sp(p,matsizes,options,orders{:});
                        
                    end
                    
                    function p=U2_U_Jz2()
                        
                        p=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        p=obj.kronall(obj.EU{3},Jz_,Jz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.z^2,obj.n.z^2
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,3,2,4],[1,3,4,2]};%[1,2,3,4], First added automatically
                        
                        p=obj.sp(p,matsizes,options,orders{:});
                        
                    end
                    
                    function p=Jz_U_U_Jzz()
                        
                        p=obj.kronall(Jz_,obj.EU{2},Jzz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z^2
                            ];
                        
                        orders={[2,1,3,4]};%[1,2,3,4], First added automatically
                        
                        p=obj.sp(p,matsizes,options,orders{:});
                        
                    end
                    
                    function p=Jz_U_U3()
                        
                        p=kron(Jz_,obj.EU{4});
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z^3,obj.n.z^3
                            ];
                        
                        orders={[2,1,3]};%[1,2,3], First added automatically
                        
                        p=obj.sp(p,matsizes,options,orders{:});
                        
                    end
                    
                    function p=Jz_U_U2_Jz()
                        
                        p=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        p=obj.kronall(Jz_,obj.EU{3},Jz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[2,1,3,4,5],...
                            [1,2,3,5,4],[2,1,3,5,4],...
                            [1,2,5,3,4],[2,1,5,3,4]};%[1,2,3,4,5], First added automatically
                        
                        p=obj.sp(p,matsizes,options,orders{:});
                        
                    end
                    
                    function p=Jz_U_U_Jz2()
                        
                        p=obj.kronall(Jz_,obj.EU{2},Jz_,Jz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[2,1,3,4,5],...
                            [1,2,4,3,5],[2,1,4,3,5],...
                            [1,2,4,5,3],[2,1,4,5,3]};%[1,2,3,4,5], First added automatically
                        
                        p=obj.sp(p,matsizes,options,orders{:});
                        
                    end
                    
                    function p=U2_Jz()
                        
                        p=kron(obj.EU{2},Jz_);
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,3,2],[3,1,2]};%[1,2,3], First added automatically
                        
                        p=obj.sp(p,matsizes,options,orders{:});
                        
                    end
                    
                end
                
                function fe=fv_a1zzzz_Jz_U2_U3_Jzz_a1zzz_U2_Jzzz()
                    
                    a1zzzz_=obj.a1zzzz{rt,rtp1};
                    
                    a1zzz_=obj.a1zzz{rt,rtp1};
                    
                    fe=obj.fv{rt,rtp1}*(a1zzzz_*kron(do_P_Jz_U2_U3(),Jzz_)*omg5+...
                        a1zzz_*kron(obj.EU{2},Jzzz_)*omg6);
                    
                    function p=do_P_Jz_U2_U3()
                        
                        options=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[2,1,3],[2,3,1]};%[1,2,3], First added automatically
                        
                        p=obj.sp(kron(Jz_,obj.EU{2}),matsizes,options,orders{:});
                        
                        if any(obj.EU{3}(:))
                            
                            p=p+obj.EU{3};
                            
                        end
                        
                    end
                    
                end
                
            end
            
            function ci=c_iterate()
                
                options=[];
                
                matsizes=[
                    obj.n.z,obj.n.z
                    obj.n.z,obj.n.z
                    obj.n.z,obj.n.z
                    obj.n.z,obj.n.z
                    obj.n.z,obj.n.z
                    ];
                
                ci=do_Jz3_EU2()+obj.kronall(Jz_,Jz_,Jz_,Jz_,Jz_)+...
                    do_Jz2_EU3()+do_Jz_EU4()+obj.EU{5};
                
                ci=W*ci*K;
                
                function p=do_Jz3_EU2()
                    
                    p=obj.kronall(Jz_,Jz_,Jz_,obj.EU{2});
                    
                    orders={[1,2,4,3,5],[1,4,2,3,5],[4,1,2,3,5],...
                        [4,1,2,5,3],[4,1,5,2,3],[4,5,1,2,3],};%[1,2,3,4,5], First added automatically
                    
                    p=obj.sp(p,matsizes,options,orders{:});
                    
                end
                
                function p=do_Jz2_EU3()
                    
                    p=sparse(0);
                    
                    if ~any(obj.EU{3}(:))
                        
                        return
                        
                    end
                    
                    p=obj.kronall(Jz_,Jz_,obj.EU{3});
                    
                    orders={[1,3,2,4,5],[1,3,4,2,5],[1,3,4,5,2],...
                        [3,1,4,5,2],[3,4,1,5,2],[3,4,5,1,2],};%[1,2,3,4,5], First added automatically
                    
                    p=obj.sp(p,matsizes,options,orders{:});
                    
                end
                
                function p=do_Jz_EU4()
                    
                    p=kron(Jz_,obj.EU{4});
                    
                    orders={[2,1,3,4,5],[2,3,1,4,5],[2,3,4,1,5],[2,3,4,5,1]};%[1,2,3,4,5], First added automatically
                    
                    p=obj.sp(p,matsizes,options,orders{:});
                    
                end
                
            end
            
        end
        
    end
    
    methods(Static)
        
    end
    
end

function C=kron_Q1_Qk_times_A(debug,varargin)

for ii=1:length(varargin)
    
    varargin{ii}=varargin{ii}.';
    
end

A=varargin{end};

varargin=varargin(1:end-1);

C=A_times_kron_Q1_Qk(debug,A,varargin{:});

C=C.';

end

function C=A_times_kron_Q1_Qk(debug,varargin)

if debug
    
    kind='direct';
    
else
    
    kind='fast';
    
end

C=utils.kronecker.A_times_kron_Q1_Qk_master(kind,varargin{:});

end

function d=store_derivatives(d,keep,iscompressed)

if iscompressed
    
    return
    
end

for ii=1:numel(d)
    
    d{ii}=d{ii}(:,keep);
    
end


end

function obj=set_compress_uncompress(obj,oo)

if ~isempty(obj.keep.v{oo})
    
    return
    
end

[obj.keep.v,obj.expand.v,obj.Comp_mat.v,obj.UnComp_mat.v]=...
    utils.kronecker.shrink_expand(obj.n.v,oo);

[obj.keep.z,obj.expand.z,obj.Comp_mat.z,obj.UnComp_mat.z]=...
    utils.kronecker.shrink_expand(obj.n.z,oo);

end

function [keep,expand,Comp_mat,UnComp_mat]=initialize_compress_expand()

proto=cell(1,5);

keep=struct(); keep.v=proto; keep.z=proto;

expand=struct(); expand.v=proto; expand.z=proto;

Comp_mat=struct(); Comp_mat.v=proto; Comp_mat.z=proto;

UnComp_mat=struct(); UnComp_mat.v=proto; UnComp_mat.z=proto;

end

function pos=set_EU_Ee_links(ne,nz,msig,order,debug)

if nargin<5
    
    debug=false;
    
end

if debug
    
    pos=proto_U_algo();
    
else
    
    pos=kron_algo();
    
    
    % pos=kron_kron_algo();
    
        
end

    function pos=kron_algo()
        
        sigpos=find(msig);
        
        col=sigpos;
        
        ii=(nz-ne+1:nz).';
        
        jj=ones(ne,1);
        
        v=jj;
        
        protorow=sparse(ii,jj,v,nz,1,ne);
        
        row=protorow;
        
        pos=cell(2,order);
        
        nek=ne;
        
        for k=1:order
            
            pos{1,k}=find(row);
            
            pos{2,k}=col*ones(nek,1);
            
            if k<order
                
                row=kron(row,protorow);
                
                col=(col-1)*nz+sigpos;
                
                nek=nek*ne;
                
            end
            
        end
        
    end

    function pos=proto_U_algo()
        
        sigpos=find(msig);
        
        ii=(nz-ne+1:nz).';
        
        jj=sigpos*ones(ne,1);
        
        protoU=sparse(ii,jj,ones(ne,1),nz,nz,ne);
        
        pos=cell(2,order);
        
        Uk=protoU;
        
        nzk=nz;
        
        for k=1:order
            
            [pos{1,k},pos{2,k}]=find(Uk);
            
            if k<order
                
                Uk=kron(Uk,protoU);
                
                nzk=nzk*nz;
                
            end
            
        end
        
        
    end

%     function pos=kron_kron_algo()
%         
%         protocol=sparse(msig);
%         
%         col=protocol;
%         
%         ii=(nz-ne+1:nz).';
%         
%         jj=ones(ne,1);
%         
%         v=jj;
%         
%         protorow=sparse(ii,jj,v,nz,1,ne);
%         
%         row=protorow;
%         
%         pos=cell(2,order);
%         
%         nek=ne;
%         
%         for k=1:order
%             
%             pos{1,k}=find(row);
%             
%             pos{2,k}=find(col)*ones(nek,1);
%             
%             if k<order
%                 
%                 row=kron(row,protorow);
%                 
%                 col=kron(col,protocol);
%                 
%                 nek=nek*ne;
%                 
%             end
%             
%         end
%         
%         
%     end

end

function [EU,maxOrder]=set_stochastic_matrices_moments(ne,nz,order,EU_Ee_links,M)

if nargin<5
    
    M=set_up_normal_moments(order,ne);
    
end

EU=cell(1,order);

nzk=nz;

failed=~true;

maxOrder=0;

for k=1:order
    
    Mk=M{k}(:);
    
    if ~failed
        
        try
            
            if any(Mk)
                
                EU{k}=sparse(EU_Ee_links{1,k},EU_Ee_links{2,k},Mk,nzk,nzk);
                
            else
                
                EU{k}=sparse(nzk,nzk);
                
            end
            
            nzk=nzk*nz;
            
            maxOrder=k;
            
        catch
            
            failed=true;
            
            warning(['order ',int2str(k),' and above fail: overly large sparse arrays'])
            
        end
        
    end
    
end

end

function M=set_up_normal_moments(order,nshocks)

% https://www.dsprelated.com/freebooks/sasp/Gaussian_Moments.html

double_factorial=@(n)prod(n:-2:0);

M=cell(1,order);

for oo=1:order
    
    M{oo}=set_up_one(oo);
    
end

    function m=set_up_one(o)
        
        if o==1
            
            m=sparse(nshocks,1);
            
        else
            
            siz=nshocks*ones(1,o);
            
            m=zeros(siz);
            
            if rem(o,2)==0
                
                e=repmat({(1:nshocks).'},1,o);
                
                IND=sub2ind(siz,e{:});
                
                m(IND)=double_factorial(o-1);
                
            end
            
        end
        
    end

end

function T=solve_generalized_sylvester(A,B,C,fplus_lambda_bf,W,gensylv_solve)

switch gensylv_solve
    
    case {1,'generalized_sylvester_direct'}
        
        T=generalized_sylvester_direct(A,B,C,fplus_lambda_bf,W);
        
    case {2,'generalized_sylvester'}
        
        T=generalized_sylvester(A,B,C,fplus_lambda_bf,W);
        
    case {3,'generalized_sylvester_direct2'}
        
        T=generalized_sylvester_direct2(A,B,C,fplus_lambda_bf,W);
        
    otherwise
        
        error(['unknown option for gensylv_solve:: ',parser.any2str(gensylv_solve)])
        
end

end

function Tzz=generalized_sylvester(A,B,C,fv_lamba_bf,W)

need_uncompress = nargin==5;

h=numel(A);

[nt,nzz]=size(A{1});

N=nt*nzz*h;

rhs=zeros(nt,nzz,h);

for ii=1:h
    
    rhs(:,:,ii)=-A{ii};
    
end

TOL=1e-6;

MAXIT=min(N,20);

solvers={@tfqmr,@bicgstab,@bicgstabl,@cgs};%,@lsqr

FLAG=1;

iter=0;

while FLAG && iter<numel(solvers)
    
    iter=iter+1;
    
    slv=solvers{iter};
    
    [Tzz0,FLAG,RELRES,ITER]  = slv(@left_hand_side,rhs(:),TOL,MAXIT);
    
    fprintf(1,'=== solver: %s, iter : %0.0f, RELRES : %0.4f === \n',func2str(slv),ITER,RELRES);
    
end

Tzz0=reshape(Tzz0,nt,nzz,h);

Tzz=cell(1,h);

for ii=1:h
    
    if need_uncompress
        
        Tzz{ii}=Tzz0(:,:,ii)*W;
        
    else
        
        Tzz{ii}=Tzz0(:,:,ii);
        
    end
    
end

    function r=left_hand_side(t)
        
        t=reshape(t,nt,nzz,h);
        
        r=zeros(nt,nzz,h);
        
        for rt=1:h
            
            r(:,:,rt)=B{rt}*t(:,:,rt);
            
            for rtp1=1:h
                
                r(:,:,rt)=r(:,:,rt)+fv_lamba_bf{rt,rtp1}*t(:,:,rtp1)*C{rt,rtp1};
                
            end
            
        end
        
        r=r(:);
        
    end

end

function Tzz=generalized_sylvester_direct(A,B,C,fv_lamba_bf,W)

need_uncompress = nargin==5;

[nt,nzz]=size(A{1});

h=numel(A);

G=zeros(nt*nzz*h);

offsetr=0;

a=zeros(nt*nzz*h,1);

for rt=1:h
    
    r_range=offsetr+(1:nt*nzz);
    
    a(r_range)=A{rt}(:);
    
    offsetc=0;
    
    for rtp1=1:h
        
        c_range=offsetc+(1:nt*nzz);
        
        G(r_range,c_range)=kron(C{rt,rtp1}.',fv_lamba_bf{rt,rtp1});
        
        if rt==rtp1
            
            G(r_range,c_range)=G(r_range,c_range)+kron(speye(nzz),B{rt});
            
        end
        
        offsetc=c_range(end);
        
    end
    
    offsetr=r_range(end);
    
end

G=sparse(G);

Tzz0=-G\a;

Tzz0=reshape(Tzz0,nt,nzz,h);

Tzz=cell(1,h);

for rt=1:h
    
    if need_uncompress
        
        Tzz{rt}=sparse(Tzz0(:,:,rt)*W);
        
    else
        
        Tzz{rt}=sparse(Tzz0(:,:,rt));
        
    end
    
end

end

function Tzz=generalized_sylvester_direct2(A,B,C,fv_lamba_bf,W)

need_uncompress = nargin==5;

[nt,nzz]=size(A{1});

h=numel(A);

offsetr=0;

a=zeros(nt*nzz*h,1);

ii=cell(h^2,1);

jj=cell(h^2,1);

ss=cell(h^2,1);

iter=0;

for rt=1:h
    
    r_range=offsetr+(1:nt*nzz);
    
    a(r_range)=A{rt}(:);
    
    offsetc=0;
    
    for rtp1=1:h
        
        iter=iter+1;
        
        c_range=offsetc+(1:nt*nzz);
        
        Grc=kron(C{rt,rtp1}.',fv_lamba_bf{rt,rtp1});
        
        if rt==rtp1
            
            Grc=Grc+kron(speye(nzz),B{rt});
            
        end
        
        % Grc(abs(Grc)<1e-9)=0;
        
        [ir,jr,ss{iter}]=find(Grc);
        
        ii{iter}=offsetr+ir;
        
        jj{iter}=offsetc+jr;
        
        offsetc=c_range(end);
        
    end
    
    offsetr=r_range(end);
    
end

ii=cell2mat(ii);

jj=cell2mat(jj);

ss=cell2mat(ss);

G=sparse(ii,jj,ss,nt*nzz*h,nt*nzz*h,numel(ss));

clear ii jj ss A C B fv_lamba_bf

Tzz0=G\(-a);

Tzz0=reshape(Tzz0,nt,nzz,h);

Tzz=cell(1,h);

for rt=1:h
    
    if need_uncompress
        
        Tzz{rt}=Tzz0(:,:,rt)*W;
        
    else
        
        Tzz{rt}=Tzz0(:,:,rt);
        
    end
    
end

end