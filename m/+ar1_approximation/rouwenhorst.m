function [Thetai,y]=rouwenhorst(mu,rho,sigma,N,p,q)
% rouwenhorst approximates and AR(1) process with a Markov chain
%
% ::
%
%   [Thetai,y]=rouwenhorst(mu,rho,sigma)
%   [Thetai,y]=rouwenhorst(mu,rho,sigma,N)
%   [Thetai,y]=rouwenhorst(mu,rho,sigma,N,p)
%   [Thetai,y]=rouwenhorst(mu,rho,sigma,N,p,q)
%
% Args:
%
%    - mu : [numeric] : unconditonal mean of the ar(1) process
%
%    - rho : [numeric] : autoregressive coefficient
%
%    - sigma : [numeric] : standard deviation of the shock
%
%    - N : [numeric|{2}] : number of states of the markov chain
%
%    - p : [numeric|{.5*(1+rho)}] : first parameter of the procedure
%
%    - q : [numeric|{.5*(1+rho)}] : second parameter of the procedure
%
% Returns:
%    :
%
%    - Thetai : [numeric] : NxN transition matrix
%
%    - y : [numeric] : vector of nodes or states
%
% Note:
%
%    - The process is assumed to be of the form
%      x_t-mu=rho*(x_{t-1}-mu)+sigma*error_term where error_term ~ N(0,1)
%
% Example:
%
%    See also:


% Reference : Karen A. Kopecky, Richard M. H. Suen (2010): "Finite state
%   Markov-chain approximations to highly persistent processes". Review of
%   Economic Dynamics 13, pp 701-714

if nargin<6
    q=.5*(1+rho);
    if nargin<5
        p=q;
        if nargin<4
            N=2;
            if nargin<3
                error('wrong number of arguments')
            end
        end
    end
end

if abs(rho)>=1
    error('process must be stationary i.e. 0<=rho<=1')
end

if p<=0||p>=1
    error('p must be such that 0<=p<=1')
end

if q<=0||q>=1
    error('q must be such that 0<=q<=1')
end

if N<=1
    error('N must greater than 1')
end

sigmaz=sigma/sqrt(1-rho^2);

psi=sqrt(N-1)*sigmaz;

y=mu+linspace(-psi,psi,N);

Theta2=[
    p,1-p
    1-q q];
Thetai=Theta2;
for ii=3:N
    Thetai=p*first_p_matrix()+...
        (1-p)*second_p_matrix()+...
        (1-q)*first_q_matrix()+...
        q*second_q_matrix();
    Thetai(2:end-1,:)=.5*Thetai(2:end-1,:);
end

    function M=first_p_matrix()
        M=zeros(ii);
        M(1:ii-1,1:ii-1)=Thetai;
    end

    function M=second_p_matrix()
        M=zeros(ii);
        M(1:ii-1,1+(1:ii-1))=Thetai;
    end

    function M=first_q_matrix()
        M=zeros(ii);
        M(1+(1:ii-1),1:ii-1)=Thetai;
    end

    function M=second_q_matrix()
        M=zeros(ii);
        M(1+(1:ii-1),1+(1:ii-1))=Thetai;
    end

end