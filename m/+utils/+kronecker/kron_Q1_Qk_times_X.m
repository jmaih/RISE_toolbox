function X = kron_Q1_Qk_times_X(X,varargin)
% X_times_kron_Q1_Qk Efficient Kronecker Multiplication of kron(Q1,Q2,...,Qk)*X
%
% ::
%
%
%   Y=kron_Q1_Qk_times_X(X,{Q1,n1},{Q2,n2},...,{Qk,nk})
%
% Args:
%
%    - **X** [matrix|vector] is a vector or a matrix whose number of rows is equal to
%      the product of the columns of all Qi's
%
%    - **varargin** [array of 2-element cells] of the form {Qi,ni}
%      - **Qi** [matrix]
%      - **ni** [integer] : number of times the kronecker of Qi appears
%
% Returns:
%    :
%
% Note:
%
%    No kronecker product is required.
%
% Example:
%
%    See also:

% Reference:
% Paulo Fernandes, Brigitte Plateau, William J. Stewart (1998): "Efficient
% Descriptor-Vector Multiplications in Stochastic Automata Networks". JACM
% 45(3): 381--414. (see page 394)
%

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