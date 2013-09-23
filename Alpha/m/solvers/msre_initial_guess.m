function T0=msre_initial_guess(A0,Aminus,Gplus01,Q,solve_initialization)
if nargin==0
    T0=struct('solve_initialization','backward');
    return
end
	n=size(A0{1},1);
	h=size(Q,1);
    T0=zeros(n,n,h);
    bkwl=any(Aminus{1},1);
    for istate=2:h
        bkwl=bkwl|any(Aminus{istate},1);
    end
    switch solve_initialization
        case {'zeros'}
        case {'backward'}
%     warnstate=warning('query','all');
%     warning('off','MATLAB:singularMatrix')
            for i_state=1:h
                % the solution that corresponds to the backward-looking model [default]
                % this should accelerate the solving and give more accurate results since
                % we strictly don't need to solve for the non-state columns i.e. the columns
                % where Aminus is zero
                T0(:,:,i_state)=-A0{i_state}\Aminus{i_state};
                if any(any(isnan(T0(:,:,i_state))))
                    T0(:,:,i_state)=-pinv(A0{i_state})*Aminus{i_state};
                    % the A0 matrix can be singular, for instance in the
                    % case of the zlb
                end
            end
%     warning(warnstate)
        case {'random'}
            level=3;
            switch level
                case 0
                    T0(:,bkwl,:)=randn(n,sum(bkwl),h);
                case 1
                    AT=zeros(n,n,h);
                    sn=solvent_norm();
                    AT(:,bkwl,:)=sn^2*randn(n,sum(bkwl),h);
                    for istate=1:h
                        QAT=zeros(n);
                        for jstate=1:h
                            QAT(:,bkwl)=QAT(:,bkwl)+Q(istate,jstate)*AT(:,bkwl,jstate);
                        end
                        QAT=QAT+A0{istate};
                        T0(:,:,istate)=-QAT\Aminus{istate};
                    end
                case 2
                    Tbkw=msre_initial_guess('back_init');
                    T0(:,bkwl,:)=Tbkw(:,bkwl,:).*rand(n,sum(bkwl),h);
                case 3
                    sn=solvent_norm();
                    T0(:,bkwl,:)=sn*randn(n,sum(bkwl),h);
            end
    end
    function sn=solvent_norm()
        n_a0=0;n_aplus=0;n_aminus=0;
        for ii_state=1:h
            n_a0=max(n_a0,norm(A0{ii_state}));
            aplus=0;
            for s1=1:h
                aplus=aplus+Gplus01{ii_state,s1};
            end
            n_aplus=max(n_aplus,norm(aplus));
            n_aminus=max(n_aminus,norm(Aminus{ii_state}));
        end
        sn=(n_a0+sqrt(n_a0^2+4*n_aplus*n_aminus))/(2*n_aplus);
    end
end
