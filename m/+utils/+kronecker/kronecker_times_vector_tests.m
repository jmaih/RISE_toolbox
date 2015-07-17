n=100;A=rand(n);B=rand(n);X=rand(n);x=X(:);
Nsim=100;
t0=zeros(2,1);
for ii=1:Nsim
    tic,z=utils.kronecker.A_kron_B_times_x(A,B,X(:),n,n);t0(1)=t0(1)+toc;
    tic,z2=kronecker_times_vector(A,B,X);t0(2)=t0(2)+toc;
end
t0=t0/min(t0);
disp(t0)