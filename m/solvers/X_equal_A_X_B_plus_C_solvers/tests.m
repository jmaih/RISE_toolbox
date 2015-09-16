solve_options=struct('TolFun',1e-7,'MaxIter',7000,'verbose',false);
%%
n=10;
T=rand(n);
% T=T*diag(rand(n,1))/T;
Q=rand(n);
Q=.5*(Q+Q');
Q=diag(diag(Q));

%%
clc
t=nan(1,4);
profile on
tic,X=sandwich_solve(T,T',Q);t(1)=toc;
profile off
profile viewer
%%
tic,[P,retcode_ds]=doubling_solve(T,T',Q,solve_options);t(2)=toc;
tic,[L,retcode_lyap]=lyapunov_equation(T,Q,solve_options);t(3)=toc;
tic,[V,retcode_tad]=sandwich_a_la_tadonki(T,T',Q,solve_options);t(4)=toc;

disp(['sandwich_solve ',num2str(max(max(abs(X-T*X*T'-Q))))])

if ~retcode_ds
    disp(['doubling_solve ',num2str(max(max(abs(P-T*P*T'-Q))))])
else
    disp(['doubling could not solve. time spent :',num2str(t(2))])
    t(2)=inf;
end

if ~retcode_lyap
    disp(['Lyapunov_solve ',num2str(max(max(abs(L-T*L*T'-Q))))])
else
    disp(['Lyapunov could not solve. time spent :',num2str(t(3))])
    t(3)=inf;
end

if ~retcode_tad
    disp(['Tadonki_solve ',num2str(max(max(abs(V-T*V*T'-Q))))])
else
    disp(['Tadonki could not solve. time spent :',num2str(t(4))])
    t(4)=inf;
end

disp(min(t)),disp(t/min(t))