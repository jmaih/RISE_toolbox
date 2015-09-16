function [V,u,stationary_vars] = schur_solve(H,GG,...
    qz_criterium,...
    lyapunov_complex_threshold)

% Solves the Lyapunov equation V =H*V*H'+ GG
%- To be implemented
% - Bartels-Stewart
% - Hessenberg-Schur

if nargin<5
    lyapunov_complex_threshold=[];
    if nargin<4
        qz_criterium = [];
    end
end
if isempty(lyapunov_complex_threshold),lyapunov_complex_threshold=1e-15;end
if isempty(qz_criterium),qz_criterium = 1+1e-6;end

nvar=size(H,1);
[V,u] = first_pass();
stationary_vars=true(1,nvar);
second_pass()

    function second_pass()
        if ~isempty(u)
            Schur_vec_tol = 1e-11;
            xx = abs(H*u);
            stationary_vars = all(xx<Schur_vec_tol,2);
            vx = nan(nvar);
            vx(stationary_vars,stationary_vars) = V(stationary_vars,stationary_vars);%aa*vx*aa'+ bb*bb'
            vx(abs(vx) < 1e-12) = 0;
            V=vx;
        end
    end

    function [V,u] = first_pass()
        [U,T] = schur(H);
        e1 = abs(ordeig(T)) > 2-qz_criterium;
        n_unstab = sum(e1);       % Number of unit roots.
        n_stab = nvar-n_unstab;  % Number of stationary variables.
        if n_unstab > 0
            % Selects stable roots
            [U,T] = ordschur(U,T,e1);
            T = T(n_unstab+1:end,n_unstab+1:end);
        end
        
        B = U(:,n_unstab+1:end)'*GG*U(:,n_unstab+1:end);
        
        Vstab = zeros(n_stab,n_stab);
        ii = n_stab;
        while ii >= 2
            if abs(T(ii,ii-1))<lyapunov_complex_threshold
                if ii == n_stab
                    f = zeros(n_stab,1);
                else
                    f = T(1:ii,:)*(Vstab(:,ii+1:end)*T(ii,ii+1:end)') + ...
                        T(ii,ii)*T(1:ii,ii+1:end)*Vstab(ii+1:end,ii);
                end
                q = eye(ii)-T(1:ii,1:ii)*T(ii,ii);
                Vstab(1:ii,ii) = q\(B(1:ii,ii)+f);
                Vstab(ii,1:ii-1) = Vstab(1:ii-1,ii)';
                ii = ii - 1;
            else
                if ii == n_stab
                    f = zeros(n_stab,1);
                    f1 = f;
                else
                    f = T(1:ii,:)*(Vstab(:,ii+1:end)*T(ii,ii+1:end)') + ...
                        T(ii,ii)*T(1:ii,ii+1:end)*Vstab(ii+1:end,ii) + ...
                        T(ii,ii-1)*T(1:ii,ii+1:end)*Vstab(ii+1:end,ii-1);
                    f1 = T(1:ii,:)*(Vstab(:,ii+1:end)*T(ii-1,ii+1:end)') + ...
                        T(ii-1,ii-1)*T(1:ii,ii+1:end)*Vstab(ii+1:end,ii-1) + ...
                        T(ii-1,ii)*T(1:ii,ii+1:end)*Vstab(ii+1:end,ii);
                end
                q = [eye(ii)-T(1:ii,1:ii)*T(ii,ii),-T(1:ii,1:ii)*T(ii,ii-1) ; ...
                    -T(1:ii,1:ii)*T(ii-1,ii),eye(ii)-T(1:ii,1:ii)*T(ii-1,ii-1) ];
                z =  q\[ B(1:ii,ii)+f ; B(1:ii,ii-1) + f1 ];
                Vstab(1:ii,ii) = z(1:ii);
                Vstab(1:ii,ii-1) = z(ii+1:end);
                Vstab(ii,1:ii-1) = Vstab(1:ii-1,ii)';
                Vstab(ii-1,1:ii-2) = Vstab(1:ii-2,ii-1)';
                ii = ii - 2;
            end
        end
        if ii == 1
            f = T(1,:)*(Vstab(:,2:end)*T(1,2:end)') + T(1,1)*T(1,2:end)*Vstab(2:end,1);
            Vstab(1,1) = (B(1,1)+f)/(1-T(1,1)*T(1,1));
        end
        V = U(:,n_unstab+1:end)*Vstab*U(:,n_unstab+1:end)';
        u = U(:,1:n_unstab);
    end

end