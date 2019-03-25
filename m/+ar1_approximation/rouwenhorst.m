%  rouwenhorst approximates and AR(1) process with a Markov chain
% 
%  ::
% 
%    [Thetai,y]=rouwenhorst(mu,rho,sigma)
%    [Thetai,y]=rouwenhorst(mu,rho,sigma,N)
%    [Thetai,y]=rouwenhorst(mu,rho,sigma,N,p)
%    [Thetai,y]=rouwenhorst(mu,rho,sigma,N,p,q)
% 
%  Args:
% 
%     mu (numeric): unconditonal mean of the ar(1) process
%     rho [numeric): autoregressive coefficient
%     sigma (numeric): standard deviation of the shock
%     N (numeric | {2}): number of states of the markov chain
%     p (numeric | {.5*(1+rho)}): first parameter of the procedure
%     q (numeric | {.5*(1+rho)}): second parameter of the procedure
% 
%  Returns:
%     :
% 
%     - **Thetai** [numeric] : NxN transition matrix
%     - **y** [numeric] : vector of nodes or states
% 
%  Note:
% 
%     - The process is assumed to be of the form
%       x_t-mu=rho*(x_{t-1}-mu)+sigma*error_term where error_term ~ N(0,1)
% 
%  Reference:
% 
%     - Karen A. Kopecky, Richard M. H. Suen (2010): "Finite state
%       Markov-chain approximations to highly persistent processes". Review of
%       Economic Dynamics 13, pp 701-714
% 
%