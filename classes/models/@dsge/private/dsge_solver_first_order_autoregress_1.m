function [Tz_pb,eigval,retcode]=dsge_solver_first_order_autoregress_1(...
    dbf_plus,ds_0,dp_0,db_0,df_0,dpb_minus,siz,options)
nregs=numel(dbf_plus);
tolerance=1e-9;
all_same=true;
if nregs>1
    ncols_bf=siz.nb+siz.nf;
    ncols_pb=siz.np+siz.nb;
    if ncols_bf
        lead_=@(x)vec(dbf_plus{x});
        lead_1=lead_(1);
    end
    curr_=@(x)vec([ds_0{x},dp_0{x},db_0{x},df_0{x}]);
    curr_1=curr_(1);
    if ncols_pb
        lag_=@(x)vec(dpb_minus{x});
        lag_1=lag_(1);
    end
    for ireg=2:nregs
        if ncols_bf
            all_same=all_same && max(abs(lead_1-lead_(ireg)))<tolerance;
        end
        all_same=all_same && max(abs(curr_1-curr_(ireg)))<tolerance;
        if ncols_pb
            all_same=all_same && max(abs(lag_1-lag_(ireg)))<tolerance;
        end
        if ~all_same
            break
        end
    end
end

if ~all_same % a diagonal transition matrix with entirely different regimes
    % solve one at a time
    retcode=0;
    nvars=siz.ns+siz.np+siz.nb+siz.nf;
    Tz_pb=nan(nvars,siz.np+siz.nb,nregs);
    eigval=cell(1,nregs);
    for ireg=1:nregs
        if ~retcode && options.occbin.do_it(ireg)
            [Tsol,eigval{ireg},retcode]=dsge_solver_first_order_autoregress_1(...
                dbf_plus(ireg),ds_0(ireg),dp_0(ireg),db_0(ireg),df_0(ireg),...
                dpb_minus(ireg),siz,options);
            if ~retcode
                Tz_pb(:,:,ireg)=Tsol;
            end
        end
    end
    if ~retcode
        eigval=cell2mat(eigval);
    end
    return
end

rise_qz_criterium=sqrt(eps);
switch lower(options.solver)
    case {'rise_1'}
        [Tz_pb,eigval,retcode]=rise_solve_constant();
    case {'klein'}
        [Tz_pb,eigval,retcode]=klein_solve();
    case {'aim'}
        [Tz_pb,eigval,retcode]=aim_solve();
    case {'sims'}
        [Tz_pb,eigval,retcode]=sims_solve();
    otherwise
        error(['unknown solver ',parser.any2str(options.solver)])
end

if nregs>1
    Tz_pb=Tz_pb(:,:,ones(1,nregs));
end

    function [Tz_pb,eigval,retcode]=aim_solve(varargin) %#ok<STOUT>
        error('the aim solver is not yet implemented')
    end

    function [Tz_pb,eigval,retcode]=sims_solve(varargin) %#ok<STOUT>
        error('the sims solver is not yet implemented')
    end

    function [TT,SS,Z,eigval,retcode]=process_eigenvalues(TT,SS,Q,Z,npred)
        % Ordered inverse eigenvalues
        %----------------------------
        eigval = ordeig(TT,SS);
        stable = abs(eigval) >= 1 + rise_qz_criterium;
        nstable = sum(stable);
        unit = abs(abs(eigval)-1) < rise_qz_criterium;
        nunit = sum(unit);
        
        retcode=0;
        if nstable+nunit<npred
            retcode=22; % no solution
        elseif nstable+nunit>npred
            retcode=21; % multiple solutions
        else
            % Clusters of unit, stable, and unstable eigenvalues.
            clusters = zeros(size(eigval));
            
            % Unit roots first.
            %------------------
            clusters(unit) = 2;
            
            % Stable roots second.
            %---------------------
            clusters(stable) = 1;
            
            % Unstable roots last.
            %---------------------
            
            % Re-order by the clusters.
            %--------------------------
            [TT,SS,~,Z] = ordqz(TT,SS,Q,Z,clusters);
        end
        % Undo the eigval inversion.
        %---------------------------
        infeigval = eigval == 0;
        eigval(~infeigval) = 1./eigval(~infeigval);
        eigval(infeigval) = Inf;
    end

    function [Tz_pb,eigval,retcode]=klein_solve()
        % put system in the form a*x(t+1)=b*x(t) where x=[x0,xf];
        %--------------------------------------------------------
        nbf=siz.nb+siz.nf;
        bf_loc=siz.ns+siz.np+(1:siz.nb+siz.nf);
        pb_loc=siz.ns+(1:siz.np+siz.nb);
        B0=[ds_0{1,1},dp_0{1,1},db_0{1,1},df_0{1,1}];
        Bminus=[zeros(siz.nd,siz.ns),dpb_minus{1,1},zeros(siz.nd,siz.nf)];
        
        a=[B0,dbf_plus{1,1}
            zeros(nbf,siz.nd+nbf)];
        a(siz.nd+1:end,bf_loc)=eye(nbf);
        
        b=[Bminus,zeros(siz.nd,nbf)
            zeros(nbf,siz.nd),-eye(nbf)];
        b=-b;
        
        [Tz_pb,eigval,retcode] = solab(siz.nd);
        if ~retcode
            Tz_pb=Tz_pb(:,pb_loc);
        end
        
        function [sol,eigval,retcode] = solab(npred)
            % npred = number of stable guys
            [TT,SS,Q,Z] = qz(full(a),full(b));      % upper triangular factorization of the matrix pencil b-za
            
            % process eigenvalues
            %---------------------
            [TT,SS,Z,eigval,retcode]=process_eigenvalues(TT,SS,Q,Z,npred);
            
            sol=[];
            if ~retcode
                z11 = Z(1:npred,1:npred);
                
                z11i = z11\eye(npred);
                s11 = TT(1:npred,1:npred);
                t11 = SS(1:npred,1:npred);
                
                dyn = s11\t11;
                sol = real(z11*dyn*z11i);
                % z21 = Z(npred+1:end,1:npred);
                % f = real(z21*z11i); % already included in the lower part of p
            end
        end
    end

    function [Tzp,eigval,retcode]=rise_solve_constant()
        % state variables (lags): pred,both
        %----------------------------------
        Apb_minus=[dpb_minus{1,1}
            sparse(siz.nb,siz.nb+siz.np)]; % auxiliary equations for 'both' variables
        
        % forward-looking variables (leads): static,both,forward
        %-------------------------------------------------------
        Asbf_plus=[sparse(siz.nd,siz.ns),dbf_plus{1,1}
            sparse(siz.nb,siz.nb+siz.nf+siz.ns)]; % auxiliary equations for 'both' variables
        
        % forward-looking variables (current): static,both,forward
        %---------------------------------------------------------
        Asbf_0=[ds_0{1,1},sparse(siz.nd,siz.nb),df_0{1,1}
            sparse(siz.nb,siz.ns),speye(siz.nb),sparse(siz.nb,siz.nf)]; % auxiliary equations for 'both' variables
        
        % state variables (current): pred,both
        %-------------------------------------
        Apb_0=[dp_0{1,1},db_0{1,1}
            sparse(siz.nb,siz.np),-speye(siz.nb)]; % auxiliary equations for 'both' variables
        [Tzp,eigval,retcode]=rise_solve_1(Asbf_plus,Apb_0,Asbf_0,Apb_minus);
        if ~retcode
            % Re-order [s,b,f,p,b] as [s,p,b,f].we simply can ignore the last b
            %------------------------------------------------------------------
            static_=1:siz.ns;
            pred_=siz.ns+siz.nb+siz.nf+(1:siz.np);
            both_=siz.ns+(1:siz.nb);
            frwrd_=siz.ns+siz.nb+(1:siz.nf);
            order_var=[static_,pred_,both_,frwrd_];
            Tzp=Tzp(order_var,:);
        end
        function [Tzp,eigval,retcode]=rise_solve_1(Afrwrd_plus,Apred_0,Afrwrd_0,Apred_minus)
            A=[Apred_0,Afrwrd_plus]; % pred,frwrd
            B=-[Apred_minus,Afrwrd_0]; % pred,frwrd
            fA=full(A); fB=full(B);
            do_norm=true;
            if do_norm
                mynorm=max([max(abs(A),[],2),max(abs(B),[],2)],[],2);
                mynorm(mynorm==0)=1;
                fA=bsxfun(@rdivide,fA,mynorm);
                fB=bsxfun(@rdivide,fB,mynorm);
            end
            npred=size(Apred_0,2);
            nfrwrd=size(Afrwrd_0,2);
            % real schur decomposition
            %-------------------------
            [TT,SS,Q,Z] = qz(fA,fB,'real');
            % so we have Q*A*Z = TT, Q*B*Z = SS.
            
            % process eigenvalues
            %---------------------
            [TT,SS,Z,eigval,retcode]=process_eigenvalues(TT,SS,Q,Z,npred);
            
            Tzp=[];
            if ~retcode
                % define
                %-------
                W=Z.';
                % partition matrices
                %-------------------
                pred=1:npred;
                frwrd=npred+(1:nfrwrd);
                W11=W(pred,pred);
                W12=W(pred,frwrd);
                W21=W(frwrd,pred);
                W22=W(frwrd,frwrd);
                
                S11=SS(pred,pred); % S12=SS(pred,frwrd); % S21=SS(frwrd,pred); % S22=SS(frwrd,frwrd);
                
                T11=TT(pred,pred); % T12=TT(pred,frwrd); % T21=TT(frwrd,pred); % T22=TT(frwrd,frwrd);
                
                % form solution: forward-looking variables
                %-----------------------------------------
                Fzp=-W22\W21;
                tmp=W11+W12*Fzp;
                
                % form solution: predetermined variables
                %---------------------------------------
                Pzp=(T11*tmp)\S11*tmp;
                
                % final solution matrix: Forward-looking+predetermined
                %-----------------------------------------------------
                Tzp=[Fzp;Pzp];
                % check again
                if any(isnan(Tzp(:)))
                    retcode=22;
                end
            end
        end
    end
end