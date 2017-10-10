function X=solve_lyapunov_equation(T,Q)

X=vartools.doubling(T,T.',Q);

end
