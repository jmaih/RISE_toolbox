function [x]=fernandez_plateau_stewart(do_left,x,varargin)
% fernandez_plateau_stewart computes kron(Q1,Q2,...,Qk)*x or x*kron(Q1,...,Qk)
%
% ::
%
%
%   Y=fernandez_plateau_stewart(do_left,X,{Q1,n1},{Q2,n2},...,{Qk,nk})
%   Y=fernandez_plateau_stewart(do_left,X,{Q1,n1},{Q2,n2},...,{Qk,nk},'-debug')
%
% Args:
%
%    - **do_left** [true|false] if true, matrix x appears to the left of the
%      kronecker product i.e. x*kron(Q1,...,Qk). If instead do_left is false,
%      matrix x appears to the right of the kronecker product i.e.
%      kron(Q1,...,Qk)*x
%
%    - **x** [matrix|vector] is a vector or a matrix whose number of columns
%      (respectively number of rows) is equal to the product of the rows
%      (respectively columns) of all Qi's
%
%    - **varargin** [array of 2-element cells|matrix] cell inputs are of the
%      form {Qi,ni}. Alternatively, an input could simply be Qi
%      - **Qi** [matrix] : must be square
%      - **ni** [integer] : number of times the kronecker of Qi appears
%
%      - **'-debug'** : entered at the end, this strings triggers the direct
%        computation of the x*kron(Q1,...,Qk) or kron(Q1,...,Qk)*x, checks the
%        differences and times each computation.
%
% Returns:
%    :
%
%    - **x** [matrix|vector] : result of the computations
%
% Note:
%
%    No kronecker product is required in the computation.
%
% Example:
%
%    x=rand(1,3^3*7^2*5^1);
%    fernandez_plateau_stewart(true,x,{rand(3),3},{rand(7),2},{rand(5),1})
%
%    See also:

% Reference:
% Paulo Fernandes, Brigitte Plateau, William J. Stewart (1998): "Efficient
% Descriptor-Vector Multiplications in Stochastic Automata Networks". JACM
% 45(3): 381--414. (see page 394)
%

debug=ischar(varargin{end});
if debug
    debug=varargin{end};
    debug(isspace(debug))=[];
    debug=strcmpi(debug,'-debug');
    if ~debug
        error('wrong debug format')
    end
    varargin=varargin(1:end-1);
    z=x;
    tic
end

nq=length(varargin);
Q=cell(1,nq);
occur=cell(1,nq);
n=cell(1,nq);
iter=0;
while ~isempty(varargin)
    iter=iter+1;
    ni=1;
    Qi=varargin{1};
    if iscell(Qi)
        ni=Qi{2};
        Qi=Qi{1};
    end
    siz=size(Qi);
    if siz(1)~=siz(2)
        error('input matrices entering the kronecker product must be square')
    end
    Q{iter}=Qi;
    n{iter}=siz(1)*ones(1,ni);
    occur{iter}=iter*ones(1,ni);
    varargin=varargin(2:end);
end
occur=cell2mat(occur);
n=cell2mat(n);
N=numel(n);

nleft=prod(n(1:N-1));
nright=1;

if do_left
    if size(x,2)~=prod(n)
        error('size of x not in synch with kronecker products')
    end
else
    if size(x,1)~=prod(n)
        error('size of x not in synch with kronecker products')
    end
end

for ii=N:-1:1
    base=0;
    jump=n(ii)*nright;
    for kk=1:nleft
        for jj=1:nright
            batch=base+jj+(0:nright:(n(ii)-1)*nright);
            if do_left
                x(:,batch)=x(:,batch)*Q{occur(ii)};
            else
                x(batch,:)=Q{occur(ii)}*x(batch,:);
            end
        end
        base=base+jump;
    end
    if ii>1
        nleft=nleft/n(ii-1);
        nright=nright*n(ii);
    end
end

Qfull=1;
if debug
    toc
    tic
    for ii=1:N
        Qfull=kron(Qfull,Q{occur(ii)});
    end
    if do_left
        ztest=z*Qfull;
    else
        ztest=Qfull*z;
    end
    max(abs(ztest(:)-x(:)))
    toc
end

end