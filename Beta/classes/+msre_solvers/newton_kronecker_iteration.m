function [T,G] = newton_kronecker_iteration(T0,Aplus,A0,Aminus,Q,nn,h,frwz)
% This algorithm expands the system X+inv(Aplus*X+A0)*Aminus=0 as in
% newton_kronecker. But instead of building a huge matrix to invert, it
% solves the problem iteratively

if nargin<8
    frwz=false;
end

% compute the G criterion
I1=eye(nn);
G0=nan(nn,nn,h);
Lst=nan(nn,nn,h);
if ~frwz
    % Create a variable that will hold Aplus*T and delete Aplus
    AplusT0=Aplus;
end
% Get the fwl and trim Aplus once and for all
% Get the bkl and trim Aminus once and for all
fwl=any(Aplus(:,:,1));
bkl=any(Aminus(:,:,1));
for istate=2:h
    fwl=fwl|any(Aplus(:,:,istate));
    bkl=bkl|any(Aminus(:,:,istate));
end
Aplus=Aplus(:,fwl,:);
Aminus=Aminus(:,bkl,:);
for st=1:h % function loop
    if frwz
        AQT=0;
        for slead=1:h
            AQT=AQT+Q(st,slead)*T0(:,:,slead);
        end
        AQT=A_times_X(Aplus(:,:,st),AQT,fwl)+A0(:,:,st); % <---         AQT=Aplus(:,:,st)*AQT+A0(:,:,st);
        
    else
        AQT=A0(:,:,st);
        for slead=1:h
            if st==1
                % this is so that those computations are done only once
                AplusT0(:,:,slead)=A_times_X(Aplus(:,:,slead),T0(:,:,slead),fwl); % <---   AplusT0(:,:,slead)=Aplus(:,:,slead)*T0(:,:,slead);
            end
            AQT=AQT+Q(st,slead)*AplusT0(:,:,slead);
        end
    end
    Lst(:,:,st)=AQT\I1;
    if any(any(isnan(Lst(:,:,st))))
        T=nan(size(T0));
        G=T;
        return
    else
        G0(:,:,st)=T0(:,:,st)+...
            X_times_A(Lst(:,:,st),Aminus(:,:,st),bkl);
    end
end
clear AplusT0 % in this case we still need Aplus and so we do not clean it.

DELTA0=zeros(nn,nn,h);
[delta,retcode,tau]=transpose_free_quasi_minimum_residual(...
    @msre_matrix_times_vector,... % coefficient matrix
    G0(:),... % right hand side
    DELTA0(:),... % initial guess
    sqrt(eps),... % tolerance level
    nn*nn*h+1,... % maximum number of iterations
    false); % display output
if retcode==201
    if tau<sqrt(eps)
        disp([mfilename,':: giving a pass without full convergence'])
        retcode=0;
    end
end

DELTA=reshape(delta,[nn,nn,h]);
G=msre_matrix_times_vector(DELTA(:));
if retcode
    T=nan(size(T0));
else
    T=T0+DELTA;
end

    function g=msre_matrix_times_vector(delta_vec)
        n2=nn^2;
        g=zeros(n2*h,1);
        for s_now=1:h
            iter_st=(s_now-1)*n2+1:s_now*n2;
            LMPRIME=transpose(...
                X_times_A(Lst(:,:,s_now),Aminus(:,:,s_now),bkl)...
                );
            if ~frwz
                g(iter_st)=-delta_vec(iter_st);
            end
            for jj=1:h
                iter_jj=(jj-1)*n2+1:jj*n2;
                if frwz
                    g(iter_st)=g(iter_st)+Q(s_now,jj)*delta_vec(iter_jj);
                else
                    g(iter_st)=g(iter_st)+Q(s_now,jj)*...
                        A_kronecker_B_times_x(LMPRIME,...
                        X_times_A(Lst(:,:,s_now),Aplus(:,:,jj),fwl),...
                        delta_vec(iter_jj),nn,nn);
%                    g(iter_st)=g(iter_st)+Q(s_now,jj)*...
%                        kronecker_times_vector(LMPRIME,...
%                        X_times_A(Lst(:,:,s_now),Aplus(:,:,jj),fwl),...
%                        reshape(delta_vec(iter_jj),nn,nn));
                end
            end
            if frwz
                g(iter_st)=-delta_vec(iter_st)+A_kronecker_B_times_x(LMPRIME,...
                    X_times_A(Lst(:,:,s_now),Aplus(:,:,s_now),fwl),...
                    g(iter_st),nn,nn);
%                g(iter_st)=-delta_vec(iter_st)+kronecker_times_vector(LMPRIME,...
%                    X_times_A(Lst(:,:,s_now),Aplus(:,:,s_now),fwl),...
%                    reshape(g(iter_st),nn,nn));
            end
        end
    end
end
