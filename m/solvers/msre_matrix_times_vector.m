function g = msre_matrix_times_vector(delta_vec,Lminus,Lplus,Q,n,h,new_version)

memory_savings=true;
% the old version is faster and solves but solves the problem
% Ap(st)EX_{t+1}+A0(st)X_{t}+Am(st)X_{t-1}+B(st)e_{t}=0. On the other hand,
% the new version is slower and solves the more realistic problem
% E[Ap(s_{t+1})X_{t+1}]+A0(st)X_{t}+Am(st)X_{t-1}+B(st)e_{t}=0. 
if nargin<7
    new_version=true;
end
n2=n^2;
g=nan(n2*h,1);

if new_version
    for st=1:h
        iter_st=(st-1)*n2+1:st*n2;
        g(iter_st)=-delta_vec(iter_st);
        if memory_savings
            LMPRIME=transpose(Lminus(:,:,1));
            Lminus=Lminus(:,:,2:end);
        else
            LMPRIME=transpose(Lminus(:,:,st));
        end
        for jj=1:h
            iter_jj=(jj-1)*n2+1:jj*n2;
            g(iter_st)=g(iter_st)+Q(st,jj)*...
                A_kronecker_B_times_x(LMPRIME,Lplus(:,:,st),delta_vec(iter_jj),n,n);
%            g(iter_st)=g(iter_st)+Q(st,jj)*...
%                kronecker_times_vector(LMPRIME,Lplus(:,:,st),reshape(delta_vec(iter_jj),n,n));
        end
    end
else
    for st=1:h
        iter_st=(st-1)*n2+1:st*n2;
        PD=0;
        for jj=1:h
            PD=PD+Q(st,jj)*delta_vec((jj-1)*n2+1:jj*n2);
        end
        g(iter_st)=-delta_vec(iter_st)+...
            A_kronecker_B_times_x(transpose(Lminus(:,:,st)),Lplus(:,:,st),PD,n,n);
%        g(iter_st)=-delta_vec(iter_st)+...
%            kronecker_times_vector(transpose(Lminus(:,:,st)),Lplus(:,:,st),reshape(PD,n,n));
    end
end


% %%
% clc
% t0=zeros(2,1);
% Nsim=500;
% n=70;
% A=rand(n);
% B=rand(n);
% X=rand(n);
% for isim=1:Nsim
%     tic,T1=kronecker_times_vector(A,B,X);t0(1)=t0(1)+toc;
%     tic,T2= kronmult({A,B},X(:));t0(2)=t0(2)+toc;
% end
% t0=t0/Nsim
% t0/min(t0)
% t0 =
% 
%     0.0007
%     0.0013
% 
% 
% ans =
% 
%     1.0000
%     2.0397