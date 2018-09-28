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
        fv_kron_Q1_Qk=@utils.kronecker.A_times_reordered_kron_Q1_Qk
        keep
        expand
        Comp_mat
        UnComp_mat
        options=struct('solve_linsyst_user_algo','',...
            'fix_point_TolFun',sqrt(eps),'fix_point_maxiter',1000,...
            'debug',false,'small_problem',true)
    end
    
    methods
        
        function obj=hops(ns,np,nb,nf,ne0,k,h,sig_in_v)
            
            if nargin
                
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
                    obj.max_order,obj.options.debug);
                
                [obj.EU,obj.max_order]=set_stochastic_matrices_moments(...
                    obj.n.e,obj.n.z,obj.max_order,obj.EU_Ee_links);
                
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
        
        function [obj,retcode]=solve(obj,f,Tz,ss,g,order,M)
                        
            if nargin<7
                
                M=[];
                
            end
            
            true_order=min(obj.max_order,order);
            
            absorb_derivatives(true_order)
            
            obj.Tz=Tz;
            
            obj.ss=ss;
            
            obj.g=g;
            
            if ~isempty(M)
                % reset the moments if the user specifies other moments
                % including e.g. skewed shocks
                obj.EU=set_stochastic_matrices_moments(obj.n.e,obj.n.z,...
                    true_order,obj.EU_Ee_links,M);
                
            end
            
            retcode=0;
            
            if true_order>1
                
                obj=set_compress_uncompress(obj,true_order);
                
                if obj.options.debug,tic,end
                
                [obj,retcode]=solve_second_order(obj);
                
                if obj.options.debug
                    
                    fprintf(1,[messenger(2,toc),'\n']);
                    
                end
                
                if true_order>2 && ~retcode
                    
                    if obj.options.debug,tic,end
                    
                    [obj,retcode]=solve_third_order(obj);
                    
                    if obj.options.debug
                        
                        fprintf(1,[messenger(3,toc),'\n']);
                        
                    end
                    
                    if true_order>3 && ~retcode
                        
                        if obj.options.debug,tic,end
                        
                        [obj,retcode]=solve_fourth_order(obj);
                        
                        if obj.options.debug
                            
                            fprintf(1,[messenger(4,toc),'\n']);
                            
                        end
                        
                        if true_order>4 && ~retcode
                            
                            if obj.options.debug,tic,end
                            
                            [obj,retcode]=solve_fifth_order(obj);
                            
                            if obj.options.debug
                                
                                fprintf(1,[messenger(5,toc),'\n']);
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
            function m=messenger(o,tt)
                
                m=sprintf('order %0.0f done in %0.4f secs, nz=%0.0f, nz^%0.0f=%0.0f',...
                    o,tt,obj.n.z,o,obj.n.z^o);
                
            end
            
            function absorb_derivatives(true_order)
                
                v=repmat('v',1,obj.max_order);
                
                for ii=1:true_order
                    
                    vi=v(1:ii);
                    
                    if ~isfield(f,vi)
                        
                        error(['derivatives of order ',int2str(ii),' not found'])
                        
                    end
                    
                    obj.(['f',vi])=f.(vi);
                    
                end
                
                clear f
                
            end
            
        end
        
    end
    
    methods(Access=private)
        
        function [obj,retcode]=solve_second_order(obj)
            
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
                    CBAR{rt,rtp1}=c_iterate_2(Jz_,obj.EU(1:2),obj.n,...
                        W,K,obj.options.small_problem,obj.options.debug);
                    
                end
                % compress
                A{rt}=A{rt}*K;
                
            end
            
            obj.B=B_; clear B_
            
            obj.fplus_lambda_bf=fplus_times_lambda_bf(obj.fv,obj.posv.bf,...
                obj.lambda.bf);
            
            [obj.Tzz,retcode]=solve_generalized_sylvester(A,obj.B,CBAR,...
                obj.fplus_lambda_bf,W,K,obj.options);
            
            function ai=a_iterate()
                
                fvv_=obj.fvv{rt,rtp1};%(:,obj.expand.v{2});
                
                ai=fvv_Evz_vz();
                
                function r=fvv_Evz_vz()
                    
                    r=A_times_kron_Q1_Qk(fvv_,a0z_,a0z_);
                    
                    r=r+A_times_kron_Q1_Qk(fvv_,a1z_,a1z_)*obj.EU{2};
                    
                end
                
            end
            
            function bi=b_iterate()
                
                fplus=obj.fv{rt,rtp1}(:,obj.posv.bf);
                
                f0=obj.fv{rt,rtp1}(:,obj.posv.spbf);
                
                Tz1=obj.Tz{rtp1}(obj.post.bf,obj.posz.pb)*obj.lambda.pb;
                
                bi=fplus*Tz1+f0;
                
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
        
        function [obj,retcode]=solve_third_order(obj)
            
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
                    
                    CBAR{rt,rtp1}=c_iterate_3(Jz_,obj.EU(1:3),obj.n,...
                        W,K,obj.options.small_problem,obj.options.debug);
                    
                end
                % compress
                A{rt}=A{rt}*K;
                
            end
            
            [obj.Tzzz,retcode]=solve_generalized_sylvester(A,obj.B,CBAR,...
                obj.fplus_lambda_bf,W,K,obj.options);
            
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
                
                TJg=utils.cr.second_order(Tzz_,Tz_,Jzz_,Jz_);
                
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
                
                ai=ai+A_times_kron_Q1_Qk(fplus_lmbd*obj.Tzz{rtp1},...
                    Jz_,Jzz_)*omg1;
                
                function fe=fvvv_Evz_vz_vz()
                    
                    fvvv_=obj.fvvv{rt,rtp1};%(:,obj.expand.v{3});
                    
                    fe=A_times_kron_Q1_Qk(fvvv_,a0z_,a0z_,a0z_);
                    
                    if any(obj.EU{3}(:))
                        
                        fe=fe+A_times_kron_Q1_Qk(fvvv_,a1z_,a1z_,...
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
                        
                        opts=[];
                        
                        fe=kron_Q1_Qk_times_A(a1z_,a1z_,obj.EU{2});
                        
                        fe=obj.fv_kron_Q1_Qk(fvvv_,matsizes,orders,opts,a0z_,fe);
                        
                    end
                    
                end
                
                function fe=fvv_Evz_vzz()
                    
                    fvv_=obj.fvv{rt,rtp1};
                    
                    EU2=obj.EU{2};
                    
                    fe=A_times_kron_Q1_Qk(fvv_,a0z_,a0zz_);
                    
                    if any(obj.EU{3}(:))
                        
                        fe=fe+A_times_kron_Q1_Qk(fvv_,a1z_,a1zz_)*obj.EU{3};
                        
                    end
                    
                    matsizes=[
                        obj.n.z,obj.n.z
                        obj.n.z,obj.n.z
                        obj.n.z,obj.n.z
                        ];
                    
                    orders={[1,3,2]};%[1,2,3], First added automatically
                    
                    opts=[];
                    
                    tmp=obj.fv_kron_Q1_Qk(A_times_kron_Q1_Qk(fvv_,a1z_,a1zz_),...
                        matsizes,orders,opts,EU2,Jz_);
                    % P(kron(EU2,Jz_))-kron(Jz_,EU2)
                                        
                    fe=fe+tmp;
                    
                    fe=fe+A_times_kron_Q1_Qk(fvv_,a0z_,a1zz_*EU2);
                    
                end
                
            end
                        
        end
        
        function [obj,retcode]=solve_fourth_order(obj)
            
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
                    
                    CBAR{rt,rtp1}=c_iterate_4(Jz_,obj.EU(1:4),obj.n,...
                        W,K,obj.options.small_problem,obj.options.debug);
                    
                end
                % compress
                A{rt}=A{rt}*K;
                
            end
            
            [obj.Tzzzz,retcode]=solve_generalized_sylvester(A,obj.B,CBAR,...
                obj.fplus_lambda_bf,W,K,obj.options);
            
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
                
                TJg=utils.cr.third_order(Tzzz_,Tzz_,Tz_,Jzzz_,Jzz_,Jz_);
                
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
                                
                fvv_=obj.fvv{rt,rtp1};
                
                ai=fvvvv_Evz_vz_vz_vz();
                
                ai=ai+fvvv_Evz2_vzz()*omg2;
                
                ai=ai+fvv_Evz_vzzz()*omg3;
                
                ai=ai+fvv_Evzz_vzz()*omg4;
                
                ai=ai+A_times_kron_Q1_Qk(fplus_lmbd*obj.Tzzz{rtp1},...
                    obj.EU{2},Jzz_)*omg2;
                
                ai=ai+A_times_kron_Q1_Qk(fplus_lmbd*obj.Tzzz{rtp1},...
                    Jz_,Jz_,Jzz_)*omg2;
                
                ai=ai+A_times_kron_Q1_Qk(fplus_lmbd*obj.Tzz{rtp1},...
                    Jzz_,Jzz_)*omg4;
                
                ai=ai+A_times_kron_Q1_Qk(fplus_lmbd*obj.Tzz{rtp1},...
                    Jz_,Jzzz_)*omg3;
                
                function fe=fvvvv_Evz_vz_vz_vz()
                    
                    fvvvv_=obj.fvvvv{rt,rtp1};%(:,obj.expand.v{4});
                    
                    EU2=obj.EU{2};
                    
                    fe=A_times_kron_Q1_Qk(fvvvv_,a0z_,a0z_,a0z_,a0z_);
                    
                    fe=fe+A_times_kron_Q1_Qk(fvvvv_,a1z_,a1z_,a1z_,a1z_)*obj.EU{4};
                    
                    opts=[];
                    
                    matsizes=[
                        obj.n.v,obj.n.z
                        obj.n.v,obj.n.z
                        obj.n.v,obj.n.z
                        obj.n.v,obj.n.z
                        ];
                    
                    fe=fe+fvvvv_do_a0z_a0z_a1z_a1z_EU2();
                    
                    fe=fe+fvvvv_do_a0z_a1zU_a1zU_a1zU();
                    
                    %----------------------------------
                    function p=fvvvv_do_a0z_a1zU_a1zU_a1zU()
                        
                        p=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        orders={[2,1,3,4],[2,3,1,4],[2,3,4,1]};%[1,2,3,4], First added automatically
                        
                        tmp=kron_Q1_Qk_times_A(a1z_,a1z_,a1z_,obj.EU{3});
                        
                        p=obj.fv_kron_Q1_Qk(fvvvv_,matsizes,orders,opts,a0z_,...
                            tmp);
                        
                    end
                    
                    function p=fvvvv_do_a0z_a0z_a1z_a1z_EU2()
                        
                        orders={[1,3,2,4],[1,3,4,2],[3,1,4,2],[3,4,1,2]};%[1,2,3,4], First added automatically
                        
                        p=kron_Q1_Qk_times_A(a1z_,a1z_,EU2);
                        
                        p=obj.fv_kron_Q1_Qk(fvvvv_,matsizes,orders,opts,...
                            a0z_,a0z_,p);
                        
                    end
                    %----------------------------------
                    
                end
                
                function fe=fvvv_Evz2_vzz()
                    
                    fvvv_=obj.fvvv{rt,rtp1};%(:,obj.expand.v{3});
                    
                    fe=A_times_kron_Q1_Qk(fvvv_,a0z_,a0z_,a0zz_);
                    
                    tmp0=kron_Q1_Qk_times_A(a1z_,a1z_,obj.EU{2});
                                                            
                    fe=fe+A_times_kron_Q1_Qk(fvvv_,tmp0,a0zz_); clear tmp0
                    
                    fvvv_I_a1zz=A_times_kron_Q1_Qk(fvvv_,speye(obj.n.v^2),a1zz_);
                    
                    fe=fe+fvvv_I_a1zz_do_a1z2_U2_U_Jz();
                    
                    fe=fe+fvvv_I_a1zz_do_a0z_a1zU_U_Jz();
                    
                    fe=fe+A_times_kron_Q1_Qk(fvvv_,a0z_,a0z_,...
                        a1zz_*obj.EU{2});
                    
                    fe=fe+fvvv_do_a0z_a1zU_a1zz_U2();
                    
                    fe=fe+A_times_kron_Q1_Qk(fvvv_,a1z_,a1z_,a1zz_)*...
                        obj.EU{4};
                    
                    function p=fvvv_I_a1zz_do_a1z2_U2_U_Jz()
                        
                        p=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        p0=kron_Q1_Qk_times_A(a1z_,a1z_,...
                            speye(obj.n.z),obj.EU{3});
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.v^2,obj.n.z^2
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,3,2]};%[1,2,3], First added automatically
                        
                        p=obj.fv_kron_Q1_Qk(fvvv_I_a1zz,matsizes,orders,opts,p0,Jz_);
                        
                    end
                    
                    function p=fvvv_I_a1zz_do_a0z_a1zU_U_Jz()
                        
                        p=kron_Q1_Qk_times_A(a1z_,speye(obj.n.z),...
                            obj.EU{2});
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[2,1,3,4],[2,1,4,3],[1,2,4,3]};%[1,2,3,4], First added automatically
                        
                        p=obj.fv_kron_Q1_Qk(fvvv_I_a1zz,matsizes,orders,...
                            opts,a0z_,p,Jz_);
                        
                    end
                    
                    function p=fvvv_do_a0z_a1zU_a1zz_U2()
                        
                        p=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z^2
                            ];
                        
                        orders={[2,1,3]};%[1,2,3], First added automatically
                        
                        tmp_=kron_Q1_Qk_times_A(a1z_,a1zz_,obj.EU{3});
                        
                        p=obj.fv_kron_Q1_Qk(fvvv_,matsizes,orders,opts,...
                            a0z_,tmp_);
                        
                    end
                    
                end
                
                function fe=fvv_Evz_vzzz()
                    
                    fe=A_times_kron_Q1_Qk(fvv_,a0z_,a1zzz_do_U2_Jz());
                    
                    if any(obj.EU{3}(:))
                        
                        fe=fe+A_times_kron_Q1_Qk(fvv_,a0z_,a1zzz_*obj.EU{3});
                        
                    end
                    
                    fe=fe+fvv_do_a1zU_a1zzz_U_Jz2();
                    
                    fe=fe+fvv_do_a1zU_a1zzz_U2_Jz();
                    
                    fe=fe+fvv_a1z_a1zzz_U4();
                    
                    fe=fe+fvv_a1z_a1zz_U2_Jzz_I_omg1();
                    
                    function p=fvv_a1z_a1zzz_U4()
                        
                        p=fvv_*kron_Q1_Qk_times_A(a1z_,a1zzz_,obj.EU{4});
                        
                    end
                    
                    function p=fvv_a1z_a1zz_U2_Jzz_I_omg1()
                        
                        p=A_times_kron_Q1_Qk(fvv_,a1z_,a1zz_);
                        
                        p=A_times_kron_Q1_Qk(p,obj.EU{2},Jzz_);
                        
                        p=A_times_kron_Q1_Qk(p,speye(obj.n.z),omg1);
                        
                    end
                    
                    function p=fvv_do_a1zU_a1zzz_U2_Jz()
                        
                        p=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,2,4,3],[1,4,2,3]};%[1,2,3,4], First added automatically
                        
                        p=A_times_kron_Q1_Qk(fvv_,a1z_,a1zzz_);
                        
                        p=obj.fv_kron_Q1_Qk(p,matsizes,orders,opts,...
                            obj.EU{3},Jz_);
                        
                    end
                    
                    function p=fvv_do_a1zU_a1zzz_U_Jz2()
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,3,2,4],[1,3,4,2]};%[1,2,3,4], First added automatically
                        
                        p0=A_times_kron_Q1_Qk(fvv_,a1z_,a1zzz_);
                        
                        p=obj.fv_kron_Q1_Qk(p0,matsizes,orders,opts,...
                            obj.EU{2},Jz_,Jz_);
                        
                    end
                    
                    function p=a1zzz_do_U2_Jz()
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,3,2],[3,1,2]};%[1,2,3], First added automatically
                        
                        p=obj.fv_kron_Q1_Qk(a1zzz_,matsizes,orders,opts,...
                            obj.EU{2},Jz_);
                    end
                    
                end
                
                function fe=fvv_Evzz_vzz()
                    
                    fe=A_times_kron_Q1_Qk(fvv_,a0zz_,a0zz_);
                    
                    fe=fe+fvv_do_a0zz_a1zzU2();
                    
                    fe=fe+fvv_do_a1zz_Jz_U_power2();
                    
                    fe=fe+fvv_do_a1zz_Jz_U_a1zz_U2();
                    
                    fe=fe+fvv_a1zz_a1zz_U4();
                    
                    function p=fvv_a1zz_a1zz_U4()
                        
                        p=A_times_kron_Q1_Qk(fvv_,a1zz_,a1zz_);
                        
                        p=p*obj.EU{4};
                        
                    end
                    
                    function p=fvv_do_a1zz_Jz_U_a1zz_U2()
                        
                        p=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        p=A_times_kron_Q1_Qk(fvv_,a1zz_,a1zz_);
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z^2,obj.n.z^2
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[3,1,2],[2,1,3],[2,3,1]};%[1,2,3], First added automatically
                        
                        p=obj.fv_kron_Q1_Qk(p,matsizes,orders,opts,...
                            Jz_,obj.EU{3});
                        
                    end
                    
                    function p=fvv_do_a1zz_Jz_U_power2()
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[2,1,3,4],[2,1,4,3],[1,2,4,3]};%[1,2,3,4], First added automatically
                        
                        p=A_times_kron_Q1_Qk(fvv_,a1zz_,a1zz_);
                        
                        p=obj.fv_kron_Q1_Qk(p,matsizes,orders,opts,...
                            Jz_,obj.EU{2},Jz_);
                        
                    end
                    
                    function p=fvv_do_a0zz_a1zzU2()
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z^2
                            obj.n.v,obj.n.z^2
                            ];
                        
                        orders={[2,1]};%[1,2], First added automatically
                        
                        p=obj.fv_kron_Q1_Qk(fvv_,matsizes,orders,opts,...
                            a0zz_,a1zz_*obj.EU{2});
                        
                    end
                    
                end
                
            end
                        
        end
        
        function [obj,retcode]=solve_fifth_order(obj)
            
            pb=obj.post.pb;
            
            bf=obj.post.bf;
            
            Jzzzz_bulk=sparse(obj.n.e*(obj.k+1)+1,obj.n.z^4);
            
            a1zzzz_bulk=sparse(obj.n.t+obj.n.pb+obj.n.e+obj.sig_in_v,obj.n.z^4);
            
            a0zzzz_bulk=sparse(obj.n.pb+obj.n.e+obj.sig_in_v,obj.n.z^4);
            
            omg1=utils.cr.omega(obj.n.z,1);
            
            omg2=utils.cr.omega(obj.n.z,2);
            
            omg3=utils.cr.omega(obj.n.z,3);
            
            % omg4=utils.cr.omega(obj.n.z,4);
            
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
                    
                    CBAR{rt,rtp1}=c_iterate_5(Jz_,obj.EU(1:5),obj.n,...
                        W,K,obj.options.small_problem,obj.options.debug);
                    
                end
                % compress
                A{rt}=A{rt}*K;
                
            end
            
            [obj.Tzzzzz,retcode]=solve_generalized_sylvester(A,obj.B,CBAR,...
                obj.fplus_lambda_bf,W,K,obj.options);
            
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
                
                TJg=utils.cr.fourth_order(Tzzzz_,Tzzz_,Tzz_,Tz_,Jzzzz_,...
                    Jzzz_,Jzz_,Jz_);
                
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
                    utils.cr.fifth_order([],Tzzzz_,Tzzz_,Tzz_,Tz_,...
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
                    
                    fvvvvv_=obj.fvvvvv{rt,rtp1};%(:,obj.expand.v{5});
                    
                    fe=A_times_kron_Q1_Qk(fvvvvv_,a0z_,a0z_,a0z_,...
                        a0z_,a0z_);
                    
                    fe=fe+fvvvvv_a0z_a1z4_U4();
                    
                    fe=fe+fvvvvv_a0z2_a1z3U3();
                    
                    fe=fe+fvvvvv_a0z3_a1z2U2();
                    
                    if any(obj.EU{5}(:))
                        
                        fe=fe+...
                            A_times_kron_Q1_Qk(fvvvvv_,a1z_,a1z_,...
                            a1z_,a1z_,a1z_)*obj.EU{5};
                        
                    end
                    
                    function fe=fvvvvv_a0z_a1z4_U4()
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v^4,obj.n.z^4
                            ];
                        
                        orders={[2,1]};%[1,2], First added automatically
                        
                        e_=kron_Q1_Qk_times_A(a1z_,a1z_,a1z_,a1z_,...
                            obj.EU{4});
                        
                        fe=obj.fv_kron_Q1_Qk(fvvvvv_,matsizes,orders,opts,...
                            a0z_,e_);
                        
                    end
                    
                    function fe=fvvvvv_a0z2_a1z3U3()
                        
                        fe=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        fe=kron_Q1_Qk_times_A(a1z_,a1z_,a1z_,obj.EU{3});
                                                
                        opts=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.v^3,obj.n.z^3
                            ];
                        
                        orders={[1,3,2],[3,1,2]};%[1,2,3], First added automatically
                        
                        fe=obj.fv_kron_Q1_Qk(fvvvvv_,matsizes,orders,opts,...
                            a0z_,a0z_,fe);
                    end
                    
                    function fe=fvvvvv_a0z3_a1z2U2()
                        
                        fe=kron_Q1_Qk_times_A(a1z_,a1z_,obj.EU{2});
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.v^2,obj.n.z^2
                            ];
                        
                        orders={[1,2,4,3],[1,4,2,3],[4,1,2,3]};%[1,2,3,4], First added automatically
                        
                        fe=obj.fv_kron_Q1_Qk(fvvvvv_,matsizes,orders,opts,...
                            a0z_,a0z_,a0z_,fe);
                        
                    end
                    
                end
                
                function fe=fvvvv_Evz3_vzz()
                    
                    fvvvv_=obj.fvvvv{rt,rtp1};%(:,obj.expand.v{4});
                    
                    fe=A_times_kron_Q1_Qk(fvvvv_,a0z_,a0z_,a0z_,...
                        a0zz_+a1zz_*obj.EU{2});
                    
                    fe=fe+fvvvv_I_a1zz_a0z2_a1zU_U_Jz_or_U2();
                    
                    fe=fe+fvvvv_a0z_a1zU_a1zU_a0zz();
                    
                    fe=fe+fvvvv_a0z_a1zU_a1zU_a1zz_U_Jz();
                    
                    fe=fe+fvvvv_a0z_a1zU_a1zU_a1zzU2();
                    
                    fe=fe+fvvvv_a1zU_a1zU_a1zU_a0zz();
                    
                    fe=fe+fvvvv_a1zU_a1zU_a1zU_a1zz_U_Jz();
                    
                    fe=fe+fvvvv_a1zU_a1zU_a1zU_a1zz_U2();
                    
                    fe=fe*omg5;
                    
                    function p=fvvvv_I_a1zz_a0z2_a1zU_U_Jz_or_U2()
                        
                        ftmp=A_times_kron_Q1_Qk(fvvvv_,speye(obj.n.v^3),a1zz_);
                        
                        p=ftmp_do_U2_Jz();
                        
                        p=p+ftmp_do_U3();
                        
                        function p=ftmp_do_U3()
                            
                            p=sparse(0);
                            
                            if ~any(obj.EU{3}(:))
                                
                                return
                                
                            end
                            
                            opts=[];
                            
                            matsizes=[
                                obj.n.v,obj.n.z
                                obj.n.v,obj.n.z
                                obj.n.v,obj.n.z
                                obj.n.z^2,obj.n.z^2
                                ];
                            
                            orders={[1,3,2,4],[3,1,2,4]};%[1,2,3,4], First added automatically
                            
                            p=kron_Q1_Qk_times_A(a1z_,...
                                speye(obj.n.z^2),obj.EU{3});
                            
                            p=obj.fv_kron_Q1_Qk(ftmp,matsizes,orders,opts,...
                            a0z_,a0z_,p);
                            
                        end
                        
                        function p=ftmp_do_U2_Jz()
                            
                            opts=[];
                            
                            matsizes=[
                                obj.n.v,obj.n.z
                                obj.n.v,obj.n.z
                                obj.n.v,obj.n.z
                                obj.n.z,obj.n.z
                                obj.n.z,obj.n.z
                                ];
                            
                            orders={[1,3,2,4,5],[3,1,2,4,5],[1,2,3,5,4],...
                                [1,3,2,5,4],[3,1,2,5,4]};%[1,2,3,4,5], First added automatically

                            p=kron_Q1_Qk_times_A(a1z_,speye(obj.n.z^2),...
                                kron(obj.EU{2},Jz_));
                            
                            p=obj.fv_kron_Q1_Qk(ftmp,matsizes,orders,opts,...
                            a0z_,a0z_,p);
                            
                        end
                        
                    end
                    
                    function p=fvvvv_a0z_a1zU_a1zU_a0zz()
                        
                        tmp=kron_Q1_Qk_times_A(a1z_,a1z_,obj.EU{2});
                                                
                        opts=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z^2
                            ];
                        
                        orders={[2,1,3,4],[2,3,1,4]};%[1,2,3], First added automatically
                        
                        p=obj.fv_kron_Q1_Qk(fvvvv_,matsizes,orders,opts,...
                            a0z_,tmp,a0zz_);
                            
                    end
                    
                    function p=fvvvv_a0z_a1zU_a1zU_a1zz_U_Jz()
                        
                        p=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        p=kron_Q1_Qk_times_A(a1z_,a1z_,eye(obj.n.z),obj.EU{3});
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[2,1,3,4,5],[2,3,1,4,5],[2,1,3,5,4],...
                            [2,3,1,5,4],[1,2,3,5,4]};%[1,2,3,4,5], First added automatically
                        
                        ftmp=A_times_kron_Q1_Qk(fvvvv_,eye(obj.n.v^3),a1zz_);
                        
                        p=obj.fv_kron_Q1_Qk(ftmp,matsizes,orders,opts,...
                            a0z_,p,Jz_);
                        
                    end
                    
                    function p=fvvvv_a0z_a1zU_a1zU_a1zzU2()
                        
                        p=kron_Q1_Qk_times_A(a1z_,a1z_,a1zz_,obj.EU{4});
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z^2
                            ];
                        
                        orders={[2,1,3,4],[2,3,1,4]};%[1,2,3,4], First added automatically
                        
                        p=obj.fv_kron_Q1_Qk(fvvvv_,matsizes,orders,opts,...
                            a0z_,p);
                        
                    end
                    
                    function p=fvvvv_a1zU_a1zU_a1zU_a0zz()
                        
                        p=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        p=kron_Q1_Qk_times_A(a1z_,a1z_,a1z_,obj.EU{3});
                        
                        p=A_times_kron_Q1_Qk(fvvvv_,p,a0zz_);
                        
                    end
                    
                    function p=fvvvv_a1zU_a1zU_a1zU_a1zz_U_Jz()
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.z^3,obj.n.z^3
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,3,2]};%[1,2,3], First added automatically
                        
                        ftmp=A_times_kron_Q1_Qk(fvvvv_,a1z_,a1z_,a1z_,a1zz_);
                        
                        p=obj.fv_kron_Q1_Qk(ftmp,matsizes,orders,opts,...
                            obj.EU{4},Jz_);
                        
                    end
                    
                    function p=fvvvv_a1zU_a1zU_a1zU_a1zz_U2()
                        
                        p=sparse(0);
                        
                        if ~any(obj.EU{5}(:))
                            
                            return
                            
                        end
                        
                        p=A_times_kron_Q1_Qk(fvvvv_,...
                            a1z_,a1z_,a1z_,a1zz_)*obj.EU{5};
                        
                    end
                    
                end
                
                function fe=fvvv_Evz_vz_vzzz()
                    
                    fvvv_=obj.fvvv{rt,rtp1};%(:,obj.expand.v{3});
                    
                    fe=A_times_kron_Q1_Qk(fvvv_,a0z_,a0z_,a0zzz_);
                    
                    fe=fe+A_times_kron_Q1_Qk(fvvv_,a0z_,a0z_,...
                        do_a1zzz_U2_Jz());
                    
                    fe=fe+fvvv_a0z_a0z_a1zzzU3();
                    
                    fe=fe+fvvv_a0z_a1zU_a1zzz_U_Jz2();
                    
                    fe=fe+fvvv_a0z_a1zU_a1zzz_U2_Jz();
                    
                    fe=fe+fvvv_a0z_a1zU_a1zzzU3();
                    
                    fe=fe+fvvv_a0z_a1zU_a1zz_U_Jzz_omg1();
                    
                    fe=fe+A_times_kron_Q1_Qk(fvvv_,...
                        kron_Q1_Qk_times_A(a1z_,a1z_,obj.EU{2}),...
                        a0zzz_);
                    
                    fe=fe+fvvv_a1zU_a1zU_a1zzz_U_Jz2();
                    
                    fe=fe+fvvv_a1zU_a1zU_a1zzz_U2_Jz();
                    
                    if any(obj.EU{5}(:))
                        
                        fe=fe+A_times_kron_Q1_Qk(fvvv_,a1z_,a1z_,a1zzz_)*obj.EU{5};
                        
                    end
                    
                    fe=fe+fvvv_a1zU_a1zU_a1zz_U_Jzz_omg1();
                    
                    fe=fe*omg6;
                    
                    function f=fvvv_a1zU_a1zU_a1zz_U_Jzz_omg1()
                        
                        f=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        f=A_times_kron_Q1_Qk(fvvv_,a1z_,a1z_,a1zz_);
                        
                        f=A_times_kron_Q1_Qk(f,obj.EU{3},Jzz_);
                        
                        f=A_times_kron_Q1_Qk(f,speye(obj.n.z^2),omg1);
                        
                    end
                    
                    function e=fvvv_a1zU_a1zU_a1zzz_U2_Jz()
                        
                        e=kron_Q1_Qk_times_A(a1z_,a1z_,speye(obj.n.z^2),obj.EU{4});
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.v^2,obj.n.z^2
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,2,4,3],[1,4,2,3]};%[1,2,3,4], First added automatically
                        
                        ftmp=A_times_kron_Q1_Qk(fvvv_,speye(obj.n.v^2),a1zzz_);
                        
                        e=obj.fv_kron_Q1_Qk(ftmp,matsizes,orders,opts,e,Jz_);
                        
                    end
                    
                    function e=fvvv_a1zU_a1zU_a1zzz_U_Jz2()
                        
                        e=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        e=kron_Q1_Qk_times_A(a1z_,a1z_,speye(obj.n.z),obj.EU{3});
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.v^2,obj.n.z^2
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,3,2,4],[1,3,4,2]};%[1,2,3,4], First added automatically
                        
                        ftmp=A_times_kron_Q1_Qk(fvvv_,...
                            speye(obj.n.v^2),a1zzz_);
                        
                        e=obj.fv_kron_Q1_Qk(ftmp,matsizes,orders,opts,...
                            e,Jz_,Jz_);
                        
                    end
                    
                    function e=do_a1zzz_U2_Jz()
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,3,2],[3,1,2]};%[1,2,3], First added automatically
                        
                        e=obj.fv_kron_Q1_Qk(a1zzz_,matsizes,orders,opts,...
                            obj.EU{2},Jz_);
                        
                    end
                    
                    function f=fvvv_a0z_a0z_a1zzzU3()
                        
                        f=sparse(0);
                        
                        if any(obj.EU{3}(:))
                            
                            f=A_times_kron_Q1_Qk(fvvv_,a0z_,a0z_,...
                                a1zzz_*obj.EU{3});
                            
                        end
                        
                    end
                    
                    function e=fvvv_a0z_a1zU_a1zzz_U_Jz2()
                        
                        e=kron_Q1_Qk_times_A(a1z_,speye(obj.n.z),obj.EU{2});
                        
                        ftmp=A_times_kron_Q1_Qk(fvvv_,speye(obj.n.v^2),a1zzz_);
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,2,4,3,5],[1,2,4,5,3],...
                            [2,1,3,4,5],[2,1,4,3,5],[2,1,4,5,3]};%[1,2,3,4,5], First added automatically
                        
                        e=obj.fv_kron_Q1_Qk(ftmp,matsizes,orders,opts,...
                            a0z_,e,Jz_,Jz_);
                        
                    end
                    
                    function e=fvvv_a0z_a1zU_a1zzz_U2_Jz()
                        
                        e=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        e=kron_Q1_Qk_times_A(a1z_,speye(obj.n.z^2),obj.EU{3});
                        
                        ftmp=A_times_kron_Q1_Qk(fvvv_,speye(obj.n.v^2),a1zzz_);
                                                
                        opts=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,2,3,5,4],[1,2,5,3,4],...
                            [2,1,3,4,5],[2,1,3,5,4],[2,1,5,3,4]};%[1,2,3,4,5], First added automatically

                        e=obj.fv_kron_Q1_Qk(ftmp,matsizes,orders,opts,...
                            a0z_,e,Jz_);
                        
                    end
                    
                    function e=fvvv_a0z_a1zU_a1zzzU3()
                        
                        e=kron_Q1_Qk_times_A(a1z_,a1zzz_,obj.EU{4});
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z^3
                            ];
                        
                        orders={[2,1,3]};%[1,2,3], First added automatically
                        
                        e=obj.fv_kron_Q1_Qk(fvvv_,matsizes,orders,opts,...
                            a0z_,e);
                        
                   end
                    
                    function e=fvvv_a0z_a1zU_a1zz_U_Jzz_omg1()
                        
                        e=kron_Q1_Qk_times_A(a1z_,speye(obj.n.z),obj.EU{2});
                        
                        ftmp=A_times_kron_Q1_Qk(fvvv_,speye(obj.n.v^2),a1zz_);
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z
                            obj.n.z^2,obj.n.z^3
                            ];
                        
                        orders={[2,1,3]};%[1,2,3], First added automatically
                        
                        e=obj.fv_kron_Q1_Qk(ftmp,matsizes,orders,opts,...
                            a0z_,e,Jzz_);
                        
                        e=A_times_kron_Q1_Qk(e,speye(obj.n.z^2),omg1);
                        
                    end
                    
                end
                
                function fe=fvvv_Evz_vzz_vzz()
                    
                    fvvv_=obj.fvvv{rt,rtp1};%(:,obj.expand.v{3});
                    
                    fe=A_times_kron_Q1_Qk(fvvv_,a0z_,a0zz_,a0zz_);
                    
                    fe=fe+fvvv_a0z_a0zz_a1zzU2();
                    
                    fe=fe+fvvv_a0z_a1zz_Jz_U_a1zzU2();
                    
                    fe=fe+fvvv_a0z_a1zz_Jz_U_a1zz_U_Jz();
                    
                    fe=fe+A_times_kron_Q1_Qk(fvvv_,a0z_,...
                        kron_Q1_Qk_times_A(a1zz_,a1zz_,obj.EU{4})...
                        );
                    
                    fe=fe+fvvv_a1zU_a1zz_U_Jz_a0zz();
                    
                    fe=fe+fvvv_a1zU_a1zzU2_a0zz();
                    
                    fe=fe+fvvv_a1zU_a1zzU2_a1zz_U_Jz();
                    
                    fe=fe+fvvv_a1zU_a1zz_U_Jz_a1zz_U_Jz();
                    
                    if any(obj.EU{5}(:))
                        
                        fe=fe+A_times_kron_Q1_Qk(fvvv_,a1z_,a1zz_,...
                            a1zz_)*obj.EU{5};
                        
                    end
                    
                    fe=fe*omg7;
                    
                    function e=fvvv_a1zU_a1zz_U_Jz_a1zz_U_Jz()
                        
                        e=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        opts=struct('skip_first',true);
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        ftmp=A_times_kron_Q1_Qk(fvvv_,a1z_,a1zz_,a1zz_);
                        
                        orders={[1,4,2,5,3],[1,4,2,3,5],[1,2,4,5,3],[1,2,4,3,5]};
                        %[1,2,3,4,5], First NOT added
                        
                        e=obj.fv_kron_Q1_Qk(ftmp,matsizes,orders,opts,...
                            obj.EU{3},Jz_,Jz_);
                        
                    end
                    
                    function e=fvvv_a1zU_a1zzU2_a1zz_U_Jz()
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.z^2,obj.n.z^2
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,2,4,3],[2,4,1,3],[1,4,2,3]};%[1,2,3,4], First added automatically
                        
                        ftmp=A_times_kron_Q1_Qk(fvvv_,a1z_,a1zz_,a1zz_);
                        
                        e=obj.fv_kron_Q1_Qk(ftmp,matsizes,orders,opts,...
                            obj.EU{4},Jz_);
                        
                    end
                    
                    function e=fvvv_a1zU_a1zzU2_a0zz()
                        
                        e=kron_Q1_Qk_times_A(a1z_,a1zz_,obj.EU{3});
                                                
                        opts=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z^2
                            obj.n.v,obj.n.z^2
                            ];
                        
                        orders={[1,3,2]};%[1,2,3], First added automatically
                        
                        e=obj.fv_kron_Q1_Qk(fvvv_,matsizes,orders,opts,...
                            e,a0zz_);
                        
                    end
                    
                    function e=fvvv_a1zU_a1zz_U_Jz_a0zz()
                        
                        opts=[];
                        
                        ftmp=A_times_kron_Q1_Qk(fvvv_,a1z_,a1zz_,...
                            speye(obj.n.v));
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.v,obj.n.z^2
                            ];
                        
                        orders={[1,3,2,4]};%[1,2,3,4], First added automatically
                        
                        e=obj.fv_kron_Q1_Qk(ftmp,matsizes,orders,opts,...
                            obj.EU{2},Jz_,a0zz_);
                        %------------
                        ftmp=A_times_kron_Q1_Qk(fvvv_,a1z_,...
                            speye(obj.n.v),a1zz_);
                        
                        % Dummy ordering
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.v,obj.n.z^2
                            obj.n.z,obj.n.z
                            ];
                        
                        % remove the dummy order
                        orders={[1,3,4,2],[1,3,2,4]};
                        
                        opts=struct('skip_first',true);
                        
                        e=e+obj.fv_kron_Q1_Qk(ftmp,matsizes,orders,opts,...
                            obj.EU{2},a0zz_,Jz_);
                        
                    end
                    
                    function e=fvvv_a0z_a1zz_Jz_U_a1zz_U_Jz()
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,2,3,5,4],[1,3,2,5,4],[1,3,2,4,5]};%[1,2,3,4,5], First added automatically
                        
                        ftmp=A_times_kron_Q1_Qk(fvvv_,...
                            speye(obj.n.v),a1zz_,a1zz_);
                        
                        e=obj.fv_kron_Q1_Qk(ftmp,matsizes,orders,opts,...
                            a0z_,Jz_,obj.EU{2},Jz_);
                        
                    end
                    
                    function e=fvvv_a0z_a1zz_Jz_U_a1zzU2()
                        
                        e=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z^2,obj.n.z^2
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,4,2,3],[1,3,2,4],[1,3,4,2]};%[1,2,3,4], First added automatically
                        
                        ftmp=A_times_kron_Q1_Qk(fvvv_,...
                            speye(obj.n.v),a1zz_,a1zz_);
                        
                        e=obj.fv_kron_Q1_Qk(ftmp,matsizes,orders,opts,...
                            a0z_,Jz_,obj.EU{3});
                        
                    end
                    
                    function e=fvvv_a0z_a0zz_a1zzU2()
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.v,obj.n.z^2
                            obj.n.v,obj.n.z^2
                            ];
                        
                        orders={[1,3,2]};%[1,2,3], First added automatically
                        
                        e=obj.fv_kron_Q1_Qk(fvvv_,matsizes,orders,opts,...
                            a0z_,a0zz_,a1zz_*obj.EU{2});
                        
                    end
                    
                end
                
                function fe=fvv_Evz_vzzzz()
                    
                    fvv_=obj.fvv{rt,rtp1};
                    
                    fe=A_times_kron_Q1_Qk(fvv_,a0z_,a0zzzz_);
                        
                    ftmp1=A_times_kron_Q1_Qk(fvv_,speye(obj.n.v),a1zzzz_);
                    
                    fe=fe+fvv_I_a1zzzz_a0z_Jz2_U2(ftmp1);
                    
                    fe=fe+fvv_I_a1zzzz_a0z_Jz_U3(ftmp1);
                    
                    fe=fe+A_times_kron_Q1_Qk(ftmp1,a0z_,obj.EU{4});
                    
                    clear ftmp1
                    
                    fe=fe+fvv_I_a1zzz_a0z_U2_Jzz_I_omg2();
                    
                    ftmp2=A_times_kron_Q1_Qk(fvv_,a1z_,a1zzzz_);
                    
                    fe=fe+fvv_a1z_a1zzzz_U_U_Jz3(ftmp2);
                    
                    fe=fe+fvv_a1z_a1zzzz_U_U2_Jz2(ftmp2);
                    
                    fe=fe+fvv_a1z_a1zzzz_U_U3_Jz(ftmp2);
                    
                    if any(obj.EU{5}(:))
                        
                        fe=fe+ftmp2*obj.EU{5};
                        
                    end
                    
                    clear ftmp2
                    
                    ftmp3=A_times_kron_Q1_Qk(fvv_,a1z_,a1zzz_);
                    
                    fe=fe+...
                        A_times_kron_Q1_Qk(...
                        fvv_a1z_a1zzz_U_U_Jz_Jzz(ftmp3)+...
                        fvv_a1z_a1zzz_U3_Jzz(ftmp3),...
                        speye(obj.n.z),omg2);
                    
                    clear ftmp3
                    
                    ftmp4=A_times_kron_Q1_Qk(fvv_,a1z_,a1zz_);
                    
                    ftmp4=A_times_kron_Q1_Qk(ftmp4,obj.EU{2},Jzzz_);
                    
                    ftmp4=A_times_kron_Q1_Qk(ftmp4,speye(obj.n.z),omg3);
                    
                    fe=fe+ftmp4;                    
                    
                    fe=fe*omg8;
                    
                    function fe=fvv_a1z_a1zzz_U3_Jzz(ftmp)
                        
                        fe=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        fe=A_times_kron_Q1_Qk(ftmp,obj.EU{3},Jzz_);
                        
                    end
                    
                    function fe=fvv_a1z_a1zzz_U_U_Jz_Jzz(ftmp)
                                                                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z^2
                            ];
                        
                        orders={[1,3,2,4]};%[1,2,3,4] automatically added
                        
                        fe=obj.fv_kron_Q1_Qk(ftmp,matsizes,orders,opts,...
                            obj.EU{2},Jz_,Jzz_);
                        
                    end
                    
                    function fe=fvv_a1z_a1zzzz_U_U3_Jz(ftmp)
                                                                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,2,3,5,4],[1,2,5,3,4],[1,5,2,3,4]};%[1,2,3,4,5] automatically added
                        
                        fe=obj.fv_kron_Q1_Qk(ftmp,matsizes,orders,opts,...
                            obj.EU{4},Jz_);
                        
                    end
                    
                    function fe=fvv_a1z_a1zzzz_U_U2_Jz2(ftmp)
                        
                        fe=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                                                                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,2,4,3,5],[1,2,4,5,3],[1,4,2,3,5],...
                            [1,4,2,5,3],[1,4,5,2,3]};%[1,2,3,4,5] automatically added
                        
                        fe=obj.fv_kron_Q1_Qk(ftmp,matsizes,orders,opts,...
                            obj.EU{3},Jz_,Jz_);
                        
                    end
                    
                    function fe=fvv_a1z_a1zzzz_U_U_Jz3(ftmp)
                                                                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,3,2,4,5],[1,3,4,2,5],[1,3,4,5,2]};%[1,2,3,4,5] automatically added
                        
                        fe=obj.fv_kron_Q1_Qk(ftmp,matsizes,orders,opts,...
                            obj.EU{2},Jz_,Jz_,Jz_);
                        
                    end
                    
                    function fe=fvv_I_a1zzz_a0z_U2_Jzz_I_omg2()
                        
                        fe=A_times_kron_Q1_Qk(fvv_,speye(obj.n.v),a1zzz_);
                        
                        fe=A_times_kron_Q1_Qk(fe,a0z_,obj.EU{2},Jzz_);
                        
                        fe=A_times_kron_Q1_Qk(fe,speye(obj.n.z),omg2);
                        
                    end
                    
                    function fe=fvv_I_a1zzzz_a0z_Jz_U3(ftmp1)
                        
                        fe=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                                                
                        opts=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,3,2,4,5],[1,3,4,2,5],[1,3,4,5,2]};%[1,2,3,4,5] automatically added
                        
                        fe=obj.fv_kron_Q1_Qk(ftmp1,matsizes,orders,opts,...
                            a0z_,Jz_,obj.EU{3});
                        
                    end
                    
                    function fe=fvv_I_a1zzzz_a0z_Jz2_U2(ftmp1)
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.v, obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,2,4,3,5],[1,2,4,5,3],[1,4,2,3,5],...
                            [1,4,2,5,3],[1,4,5,2,3]};%[1,2,3,4,5] automatically added
                        
                        fe=obj.fv_kron_Q1_Qk(ftmp1,matsizes,orders,opts,...
                            a0z_,Jz_,Jz_,obj.EU{2});
                        
                    end
                    
                end
                
                function fe=fvv_Evzz_vzzz()
                    
                    fvv_=obj.fvv{rt,rtp1};
                    
                    fe=A_times_kron_Q1_Qk(fvv_,a0zz_,a0zzz_);
                    
                    fe=fe+A_times_kron_Q1_Qk(fvv_,a1zz_*obj.EU{2},a0zzz_);
                    
                    ftmp1=A_times_kron_Q1_Qk(fvv_,a1zz_,a1zzz_);
                    
                    fe=fe+fvv_a1zz_a1zzz_Jz_U_U_Jz2(ftmp1);
                    
                    fe=fe+fvv_a1zz_a1zzz_U2_U_Jz2(ftmp1);
                    
                    fe=fe+fvv_I_a1zzz_a0zz_Jz_U2();
                    
                    ftmp2=A_times_kron_Q1_Qk(fvv_,a1zz_,a1zzz_);
                    
                    fe=fe+fvv_a1zz_a1zzz_Jz_U_U2_Jz(ftmp2);
                    
                    fe=fe+fvv_a1zz_a1zzz_U2_U2_Jz(ftmp2); clear ftmp2
                    
                    if any(obj.EU{3}(:))
                        
                        fe=fe+A_times_kron_Q1_Qk(fvv_,a0zz_,...
                            a1zzz_*obj.EU{3});
                        
                    end
                    
                    fe=fe+fvv_a1zz_a1zzz_Jz_U_U3(ftmp1);
                    
                    if any(obj.EU{5}(:))
                        
                        fe=fe+ftmp1*obj.EU{5};
                        
                    end, clear ftmp1
                    
                    fe=fe+fvv_a1zz_a1zz_Jz_U_U_Jzz_I_omg1();
                    
                    if any(obj.EU{3}(:))
                        
                        ftmp=A_times_kron_Q1_Qk(fvv_,a1zz_,a1zz_);
                        
                        ftmp=A_times_kron_Q1_Qk(ftmp,obj.EU{3},Jzz_);
                        
                        ftmp=A_times_kron_Q1_Qk(ftmp,speye(obj.n.z^2),omg1);
                        
                        fe=fe+ftmp; clear ftmp
                        
                    end
                                        
                    fe=fe*omg9;
                    
                    function fe=fvv_a1zz_a1zz_Jz_U_U_Jzz_I_omg1()
                        
                        fe=A_times_kron_Q1_Qk(fvv_,a1zz_,a1zz_);
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z^2
                            ];
                        
                        orders={[2,1,3,4]};%[1,2,3,4], First added automatically
                        
                        fe=obj.fv_kron_Q1_Qk(fe,matsizes,orders,opts,...
                            Jz_,obj.EU{2},Jzz_);
                        
                        fe=A_times_kron_Q1_Qk(fe,speye(obj.n.z^2),omg1);
                        
                    end
                    
                    function fe=fvv_a1zz_a1zzz_Jz_U_U3(ftmp)
                                                
                        opts=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z^3,obj.n.z^3
                            ];
                        
                        orders={[2,1,3]};%[1,2,3], First added automatically
                        
                        fe=obj.fv_kron_Q1_Qk(ftmp,matsizes,orders,opts,...
                            Jz_,obj.EU{4});
                        
                    end
                    
                    function fe=fvv_a1zz_a1zzz_U2_U2_Jz(ftmp)
                                                
                        opts=[];
                        
                        matsizes=[
                            obj.n.z^2,obj.n.z^2
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,2,4,3],[1,4,2,3]};%[1,2,3,4], First added automatically
                        
                        fe=obj.fv_kron_Q1_Qk(ftmp,matsizes,orders,opts,...
                            obj.EU{4},Jz_);
                        
                    end
                    
                    function fe=fvv_a1zz_a1zzz_Jz_U_U2_Jz(ftmp)
                        
                        fe=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                                                
                        opts=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,2,3,5,4],[1,2,5,3,4],[2,1,3,4,5],...
                            [2,1,3,5,4],[2,1,5,3,4]};%[1,2,3,4,5], First added automatically
                        
                        fe=obj.fv_kron_Q1_Qk(ftmp,matsizes,orders,opts,...
                            Jz_,obj.EU{3},Jz_);
                        
                    end
                    
                    function fe=fvv_I_a1zzz_a0zz_Jz_U2()
                        
                        fe=A_times_kron_Q1_Qk(fvv_,...
                            speye(obj.n.v),a1zzz_);
                                                
                        opts=[];
                        
                        matsizes=[
                            obj.n.v,obj.n.z^2
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,3,2,4],[1,3,4,2]};%[1,2,3,4], First added automatically
                        
                        fe=obj.fv_kron_Q1_Qk(fe,matsizes,orders,opts,...
                            a0zz_,Jz_,obj.EU{2});
                        
                    end
                    
                    function fe=fvv_a1zz_a1zzz_U2_U_Jz2(ftmp)
                        
                        fe=sparse(0);
                        
                        if ~any(obj.EU{3}(:))
                            
                            return
                            
                        end
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.z^2,obj.n.z^2
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,3,2,4],[1,3,4,2]};%[1,2,3,4], First added automatically
                        
                        fe=obj.fv_kron_Q1_Qk(ftmp,matsizes,orders,opts,...
                            obj.EU{3},Jz_,Jz_);
                        
                    end
                    
                    function fe=fvv_a1zz_a1zzz_Jz_U_U_Jz2(ftmp)
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            ];
                        
                        orders={[1,2,4,3,5],[1,2,4,5,3],[2,1,3,4,5],...
                            [2,1,4,3,5],[2,1,4,5,3]};%[1,2,3,4,5], First added automatically
                        
                        fe=obj.fv_kron_Q1_Qk(ftmp,matsizes,orders,opts,...
                            Jz_,obj.EU{2},Jz_,Jz_);
                        
                    end
                    
                end
                
                function fe=fv_a1zzzz_Jz_U2_U3_Jzz_a1zzz_U2_Jzzz()
                    
                    a1zzzz_=obj.a1zzzz{rt,rtp1};
                    
                    a1zzz_=obj.a1zzz{rt,rtp1};
                    
                    fv_=obj.fv{rt,rtp1};
                    
                    ftmp=fv_*a1zzzz_;
                    
                    fe1=fv_a1zzzz_Jz_U2_Jzz(ftmp);
                    
                    if any(obj.EU{3}(:))
                        
                        fe1=fe1+A_times_kron_Q1_Qk(ftmp,obj.EU{3},Jzz_);
                        
                    end
                    
                    clear ftmp
                    
                    fe2=A_times_kron_Q1_Qk(fv_*a1zzz_,obj.EU{2},Jzzz_);
                                                            
                    fe=fe1*omg5+fe2*omg6;
                    
                    function fe=fv_a1zzzz_Jz_U2_Jzz(ftmp)
                        
                        opts=[];
                        
                        matsizes=[
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z
                            obj.n.z,obj.n.z^2
                            ];
                        
                        orders={[2,1,3,4],[2,3,1,4]};%[1,2,3], First added automatically
                                                 
                        fe=obj.fv_kron_Q1_Qk(ftmp,matsizes,orders,opts,...
                            Jz_,obj.EU{2},Jzz_);
                        
                   end
                                        
                end
                
            end
            
        end
        
    end
    
    methods(Static,Hidden)
        
        varargout=set_up_normal_moments(varargin)
        
        varargout=set_up_moments(varargin)
        
    end
    
end

function ci=c_iterate_2(Jz_,EU,n,W,K,store,debug)

if store
    
    ci=engine(W,K);
    
else
    
    clear W K
    
    ci=@engine;
    
end

    function ci=engine(W,K)
        
        ci=A_times_kron_Q1_Qk(W,Jz_,Jz_);
        
        ci=ci+W*EU{2};
        
        ci=ci*K;
        
    end

end

function ci=c_iterate_3(Jz_,EU,n,W,K,store,debug)

if store
    
    ci=engine(W,K);
    
else
    
    clear W K
    
    ci=@engine;
    
end

    function WciK=engine(W,K)
        
        f_kron_J_U=@utils.kronecker.A_times_reordered_kron_Q1_Qk;
        
        matsizes=[
            n.z,n.z
            n.z,n.z
            n.z,n.z
            ];
        
        orders={[2,1,3],[2,3,1]};%[1,2,3], First added automatically
        
        opts=[];
        
        WciK=f_kron_J_U(W,matsizes,orders,opts,Jz_,EU{2});
        
        WciK=WciK+W*EU{3};
        
        WciK=WciK+A_times_kron_Q1_Qk(W,Jz_,Jz_,Jz_);
        
        WciK=WciK*K;
        
    end

end

function ci=c_iterate_4(Jz_,EU,n,W,K,store,debug)

if store
    
    ci=engine(W,K);
    
else
    
    clear W K
    
    ci=@engine;
    
end

    function ci=engine(W,K)
        
        f_kron_J_U=@utils.kronecker.A_times_reordered_kron_Q1_Qk;
        
        opts=[];
        
        matsizes=[
            n.z,n.z
            n.z,n.z
            n.z,n.z
            n.z,n.z
            ];
        
        ci=A_times_kron_Q1_Qk(W,Jz_,Jz_,Jz_,Jz_);
        
        ci=ci+W_do_Jz_Jz_EU2();
        
        ci=ci+W_do_Jz_EU3();
        
        ci=ci+W*EU{4};
        
        ci=ci*K;
        
        function p=W_do_Jz_Jz_EU2()
            
            orders={[1,3,2,4],[1,3,4,2],[3,1,4,2],[3,4,1,2]};%[1,2,3,4], First added automatically
            
            p=f_kron_J_U(W,matsizes,orders,opts,...
                Jz_,Jz_,EU{2});
            
        end
        
        function p=W_do_Jz_EU3()
            
            p=sparse(0);
            
            if ~any(EU{3}(:))
                
                return
                
            end
            
            orders={[2,1,3,4],[2,3,1,4],[2,3,4,1]};%[1,2,3,4], First added automatically
            
            p=f_kron_J_U(W,matsizes,orders,opts,...
                Jz_,EU{3});
            
        end
        
    end

end
            
function ci=c_iterate_5(Jz_,EU,n,W,K,store,debug)

if store
    
    ci=engine(W,K);
    
else
    
    clear W K
    
    ci=@engine;
    
end

    function ci=engine(W,K)
        
        f_kron_J_U=@utils.kronecker.A_times_reordered_kron_Q1_Qk;
        
        opts=[];
        
        matsizes=[
            n.z,n.z
            n.z,n.z
            n.z,n.z
            n.z,n.z
            n.z,n.z
            ];
        
        ci=W_Jz3_U2();
        
        ci=ci+A_times_kron_Q1_Qk(W,Jz_,Jz_,Jz_,Jz_,Jz_);
        
        ci=ci+W_Jz2_U3();
        
        ci=ci+W_Jz_U4();
        
        ci=ci+W*EU{5};
        
        ci=ci*K;
        
        function c=W_Jz_U4()
            
            if debug
                
                orders=set_orders({'J',1},{'U',4});
                
                % the first order is added automatically
                orders=orders(2:end);
                
            else
                
                orders={[2,1,3,4,5],[2,3,1,4,5],[2,3,4,1,5],...
                    [2,3,4,5,1]};
                
            end
            
            c=f_kron_J_U(W,matsizes,orders,opts,Jz_,EU{4});
            
        end
        
        function c=W_Jz2_U3()
            
            c=sparse(0);
            
            if ~any(EU{3}(:))
                
                return
                
            end
            
            if debug
                
                orders=set_orders({'J',2},{'U',3});
                
                % the first order is added automatically
                orders=orders(2:end);
                
            else
                
                orders={[1,3,2,4,5],[1,3,4,2,5],[1,3,4,5,2],...
                    [3,1,2,4,5],[3,1,4,2,5],[3,1,4,5,2],...
                    [3,4,1,2,5],[3,4,1,5,2],[3,4,5,1,2]};
                
            end
            
            c=f_kron_J_U(W,matsizes,orders,opts,...
                Jz_,Jz_,Jz_,EU{2});
            
        end
        
        function c=W_Jz3_U2()
            
            if debug
                
                orders=set_orders({'J',3},{'U',2});
                
                % the first order is added automatically
                orders=orders(2:end);
                
            else
                
                orders={[1,2,4,3,5],[1,2,4,5,3],[1,4,2,3,5],...
                    [1,4,2,5,3],[1,4,5,2,3],[4,1,2,3,5],...
                    [4,1,2,5,3],[4,1,5,2,3],[4,5,1,2,3]};
                
            end
            
            c=f_kron_J_U(W,matsizes,orders,opts,...
                Jz_,Jz_,Jz_,EU{2});
            
        end
        
    end

end

function o=set_orders(varargin)

[~,~,o]=utils.gridfuncs.find_combosx(varargin{:});

[r,c]=size(o);

o=mat2cell(o,ones(r,1),c);

end

function C=kron_Q1_Qk_times_A(varargin)

for ii=1:length(varargin)
    
    varargin{ii}=varargin{ii}.';
    
end

A=varargin{end};

varargin=varargin(1:end-1);

C=A_times_kron_Q1_Qk(A,varargin{:});

C=C.';

end

function C=A_times_kron_Q1_Qk(varargin)

C=utils.kronecker.A_times_kron_Q1_Qk_master([],varargin{:});

end

function obj=set_compress_uncompress(obj,oo)

if ~isempty(obj.keep.z{oo})
    
    return
    
end

[obj.keep.z,obj.expand.z,obj.Comp_mat.z,obj.UnComp_mat.z]=...
    utils.kronecker.shrink_expand(obj.n.z,oo);

end

function [keep,expand,Comp_mat,UnComp_mat]=initialize_compress_expand()

proto=cell(1,5);

keep=struct(); keep.z=proto;

expand=struct(); expand.z=proto;

Comp_mat=struct(); Comp_mat.z=proto;

UnComp_mat=struct(); UnComp_mat.z=proto;

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

function [EU,maxOrder]=set_stochastic_matrices_moments(ne,nz,order,...
    EU_Ee_links,M)

if nargin<5
    
    M=[];
    
end

if isempty(M)

    M=hops.set_up_normal_moments(order,ne);
    
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
                
                % trim the fat: the Mk matrices are very very sparse
                %----------------------------------------------------
                good=find(Mk);
                
                EU{k}=sparse(EU_Ee_links{1,k}(good),...
                    EU_Ee_links{2,k}(good),Mk(good),nzk,nzk);
                
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

function [T,retcode]=solve_generalized_sylvester(A,B,C,fplus_lambda_bf,W,K,opts)

% This seems to be crucial
%-------------------------
precondition()

[T,retcode]=generalized_sylvester(A,B,C,fplus_lambda_bf,W,K,opts);
        
    function precondition()
        
        nt=size(A{1},1);
        
        I=speye(nt);
        
        h=numel(B);
        
        for ii=1:h
            
            iB=B{ii}\I;
            
            A{ii}=iB*A{ii};
            
            for jj=1:h
                
                fplus_lambda_bf{ii,jj}=iB*fplus_lambda_bf{ii,jj};
                
            end
            
            B{ii}=I;
            
        end
        
    end

end

function [Tzz,retcode]=generalized_sylvester(A,B,C,fv_lamba_bf,W,K,opts)

func_type=isa(C{1,1},'function_handle');

need_uncompress = ~isempty(W);

h=numel(A);

[nt,nzz]=size(A{1});

rhs=zeros(nt,nzz,h);

for ii=1:h
    
    rhs(:,:,ii)=-A{ii};
    
end

Tzz0=[];

[Tzz0,retcode]=utils.optim.linear_systems_solver(@left_hand_side,rhs(:),...
                Tzz0,opts);

Tzz0=reshape(Tzz0,nt,nzz,h);

Tzz=cell(1,h);

for ii=1:h
    
    if need_uncompress
        
        Tzz{ii}=sparse(Tzz0(:,:,ii)*W);
        
    else
        
        Tzz{ii}=sparse(Tzz0(:,:,ii));
        
    end
    
end

    function r=left_hand_side(t)
        
        t=reshape(t,nt,nzz,h);
        
        r=zeros(nt,nzz,h);
        
        for rt=1:h
            
            r(:,:,rt)=B{rt}*t(:,:,rt);
            
            for rtp1=1:h
                
                FT=sparse(fv_lamba_bf{rt,rtp1}*t(:,:,rtp1));
                
                if func_type
                    
                    r(:,:,rt)=r(:,:,rt)+C{rt,rtp1}(FT*W,K);
                    
                else
                    
                    r(:,:,rt)=r(:,:,rt)+FT*C{rt,rtp1};
                    
                end
                
            end
            
        end
        
        r=r(:);
        
    end

end