function [Aplus,A0,Aminus,B,Q,T,R]=generate_markov_switching_rational_expectations_problem(n,h,x)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

if nargin<3
    x=0;
    if nargin<2
        h=1;
        if nargin<1
            n=1;
        end
    end
end
Q=rand(h);
Q=bsxfun(@rdivide,Q,sum(Q,2));
Q(:,end)=1-sum(Q(:,end-1),2);
T=rand(n,n,h);
Aplus=rand(n,n,h);
A0=rand(n,n,h);
B=rand(n,x,h);
R=rand(n,x,h);
Aminus=nan(n,n,h);
for st=1:h
    ATA=0;
    for j=1:h
        ATA=ATA+Q(st,j)*T(:,:,j);
    end
    ATA=Aplus(:,:,st)*ATA+A0(:,:,st);
    Aminus(:,:,st)=-ATA*T(:,:,st);
    R(:,:,st)=-ATA\B(:,:,st);
end