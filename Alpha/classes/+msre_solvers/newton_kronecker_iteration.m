function [T,G] = newton_kronecker_iteration(T0,Gplus01,A0,Aminus,~,nn,h,~)
% This algorithm expands the system X+inv(Aplus*X+A0)*Aminus=0 as in
% newton_kronecker. But instead of building a huge matrix to invert, it
% solves the problem iteratively

% compute the G criterion
I1=eye(nn);
G0=nan(nn,nn,h);
Lst=cell(1,h);

for st=1:h % function loop
    AQT=A0{st};
    for slead=1:h
        % AQT=AQT+Q(st,slead)*Aplus{st}*T0(:,:,slead);
        AQT=AQT+Gplus01{st,slead}*T0(:,:,slead);
    end
    Lst{st}=AQT\I1;
    if any(isnan(Lst{st}(:)))
        T=nan(size(T0));
        G=T;
        return
    else
        G0(:,:,st)=T0(:,:,st)+Lst{st}*Aminus{st};
    end
end

delta0=zeros(nn*nn*h,1);
[delta,retcode,tau]=transpose_free_quasi_minimum_residual(...
    @msre_matrix_times_vector,... % coefficient matrix
    G0(:),... % right hand side
    delta0,... % initial guess
    sqrt(eps),... % tolerance level
    nn*nn*h+1,... % maximum number of iterations
    false); % display output
if retcode==201
    if tau<sqrt(eps)
        disp([mfilename,':: giving a pass without full convergence'])
        retcode=0;
    end
end

G=msre_matrix_times_vector(delta);
if retcode
    T=nan(size(T0));
else
	DELTA=reshape(delta,[nn,nn,h]);
    T=T0+DELTA;
end

    function g=msre_matrix_times_vector(delta_vec)
        n2=nn^2;
        g=zeros(n2*h,1);
        for s0=1:h
            iter_s0=(s0-1)*n2+1:s0*n2;
            LMPRIME=transpose(Lst{s0}*Aminus{s0});
            for s1=1:h
                iter_s1=(s1-1)*n2+1:s1*n2;
				Ls0_s1=Lst{s0}*Gplus01{s0,s1};
%                g(iter_s0)=g(iter_s0)+kron(LMPRIME,Ls0_s1)*delta_vec(iter_s1);
                g(iter_s0)=g(iter_s0)+A_kronecker_B_times_x(LMPRIME,Ls0_s1,delta_vec(iter_s1),nn,nn);
            end
            g(iter_s0)=g(iter_s0)-delta_vec(iter_s0);
        end
    end
end
