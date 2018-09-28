clc,clear,close all

p=5;

lambda=rand;

max_size=10;

n=zeros(1,p);

A=cell(1,p);

for ii=1:p
    
    ni=0;
    
    while ni<=3||(ii>1 && any(ni==n(1:ii-1)))
        
        ni=randi(max_size,1);
        
    end
    
    n(ii)=ni;
    
    A{ii}=randn(n(ii));
        
end

N=prod(n);

x=randn(N,1);

b=utils.kronecker.kron_Q1_Qk_times_A(x,A{:})-lambda*x; % <-- b=(kA-lambda*speye(N))*x;

%%
clc

tic,x1=utils.kronecker.kp_shift_solve(A,lambda,b,false);toc

tic,x2=utils.kronecker.kp_shift_solve(A,lambda,b,true);toc

N

max([max(abs(x1-x)),max(abs(x2-x))])