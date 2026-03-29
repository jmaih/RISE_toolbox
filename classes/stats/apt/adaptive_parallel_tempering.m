%  Also known as:
%  - Replica exchange Monte Carlo
%  - Metropolis coupled Markov chain Monte Carlo
% 
%  Out:
%    X      -- Simulated variables in 3D array (dim*levels*N)
%    m      -- The final proposal means (vectors of length dim)
%    Rchol      -- The final proposal covariances (dim*dim matrices)
%    log_c  -- The final scaling factors (vector of length H)
%    lambda   -- The final inverse temperatures (vector of length H)
%    stats  -- Acceptance rate statistics etc.
%