function X=solve_lyapunov_equation(T,Q)
% INTERNAL FUNCTION
%

X=vartools.doubling(T,T.',Q);

end
