var R, DELTA, Z, C, W, N, K, J, PSTAR, PI, MC;
varexo EZ, ER;

parameters Beta, sigma, psi, epsilon, theta, kappa_c, kappa_pi, rho_r, rho, PI_ss, DELTA_ss, Z_ss, N_ss, C_ss,
    MC_ss, W_ss, phi, K_ss, J_ss, PSTAR_ss, R_ss;

Beta = 0.99;
sigma = 2;
psi = 3;
epsilon = 6;
theta = 0.75;
kappa_c = 0.5;
kappa_pi = 1.5;
rho_r = 0.8;
rho = 0.8;
PI_ss = 1;
DELTA_ss = 1;
Z_ss = 1;
N_ss = 1/3;
C_ss = Z_ss*N_ss;
MC_ss = (epsilon-1)/epsilon;
W_ss = C_ss*MC_ss/N_ss;
phi = W_ss/(N_ss^psi*C_ss^sigma);
K_ss = C_ss^(1-sigma)*MC_ss/(1 - theta*Beta*PI_ss^epsilon);
J_ss = C_ss^(1-sigma)/(1 - theta*Beta*PI_ss^(epsilon-1));
PSTAR_ss = (epsilon/(epsilon - 1))*(K_ss/J_ss);
R_ss = PI_ss/Beta;

model;

    Beta*R*((C/C(1))^sigma*(1/PI(1))) - 1;
    W - phi*N^psi*C^sigma;
    C - Z*N/DELTA;
    PSTAR - (epsilon/(epsilon - 1))*(K/J);
    K - C^(1-sigma)*MC - theta*Beta*(PI(1))^epsilon*K(1);
    J - C^(1-sigma) - theta*Beta*PI(1)^(epsilon-1)*J(1);
    1 - (1 - theta)*PSTAR^(1-epsilon) - theta*PI^(epsilon-1);
    DELTA - (1-theta)*PSTAR^(-epsilon) - theta*(PI)^epsilon*DELTA(-1);
    log(R/R_ss) - (1-rho_r)*(kappa_c*log(C/C_ss) + kappa_pi*log(PI/PI_ss)) - rho_r*log(R(-1)/R_ss) - ER;
    MC - W*N/C;
    log(Z/Z_ss) - rho*log(Z(-1)/Z_ss) - EZ;


end;

initval;

R = R_ss;
DELTA = DELTA_ss;
Z = Z_ss;
C = C_ss;
W = W_ss;
N = N_ss;
K = K_ss;
J = J_ss;
PSTAR = PSTAR_ss;
PI = PI_ss;
MC = MC_ss;

end;

shocks;

var EZ; 
stderr 0.01;
var ER;
stderr 0.01;

end;

%stoch_simul(irf = 40, order = 1);
stoch_simul(irf = 40, order = 3, k_order_solver);

















