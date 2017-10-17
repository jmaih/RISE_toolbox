function [S,myeig,retcode]=rise_cyclic_reduction(Aplus,A0,Aminus,Q,T0,tol,maxiter)

if size(Q,1)>1||size(T0,3)>1
    
    error('Cyclic reduction solver not appropriate for regime switching models')
    
end

myeig=[];

debug=[];

[S,retcode]=msre_solvers.cyclic_reduction(Aplus,A0,Aminus,tol,maxiter,debug);

end