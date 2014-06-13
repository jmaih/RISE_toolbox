function X = X_times_kron_Q1_Qk(X,varargin)
% X_TIMES_KRON_Q1_Qk Efficient Kronecker Multiplication
%
% Y=X_times_kron_Q1_Qk(X,{Q1,n1},{Q2,n2},...,{Qk,nk}) computes
%     Y = (kron^n1(Q1) kron kron^n2(Q2) ...  kron kron^nk(Qk))*X
% without ever using a single kronecker.
% The inputs are such that :
%                          Qi is a square matrix and occurs ni times.
%                          X is a vector or a matrix whose number of rows
%                          is equal to the product of the rows of all Qi's
% This code implements the algorithm from page 394 of Paulo Fernandes,
% Brigitte Plateau, William J. Stewart (1998): "Efficient Descriptor-Vector
% Multiplications in Stochastic Automata Networks". JACM 45(3): 381--414
% (doi:10.1145/278298.278303). 
%
% See also KRON
%
% Example:
% n=[10,40];
% calls=[3,1];
% nx=prod(n.^calls);
% x = randn(nx,1); % generate data
% A1={randn(n(1)),calls(1)};
% A2={randn(n(2)),calls(2)};
% tic,y = X_times_kron_Q1_Qk(x,A1,A2);toc
% tic
% Qfull = 1;
% for ii=1:calls(1)
%     Qfull=kron(Qfull,A1{1});
% end
% for ii=1:calls(2)
%     Qfull=kron(Qfull,A2{1});
% end
% toc
% z = Qfull*x;
% max(abs(y-z))

nmat=length(varargin);
occur=cell(1,nmat);
for imat=1:nmat
    occur{imat}=imat*ones(1,varargin{imat}{2});
end
occur=cell2mat(occur);

N = length(occur);    % number of matrices
n = zeros(N,1);
nright = 1;
for ii=1:N-1
    n(ii) = size(varargin{occur(ii)}{1},1);
end
n(N) = size(varargin{occur(N)}{1},1);
nleft =prod(n(1:N-1));

for ii=N:-1:1
    base = 0;
    jump = n(ii)*nright;
    for k=1:nleft
        for jj=1:nright
            index = base+jj+(0:nright:nright*(n(ii)-1));
            X(index,:) = varargin{occur(ii)}{1}*X(index,:);
        end
        base = base+jump;
    end
    nleft = nleft/n(max(ii-1,1));
    nright = nright*n(ii);
end