function Y = kron_times_vector(T1,T2,X)
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

% this function computes kron(T1,T2)*vec(X)
% kron(T1,T2)*U = [kron(T1,I)*kron(I,T2)]*U=kron(T1,I)*[kron(I,T2)*U]
W = kron_it(T2,X);
Y = kron_ti(T1,W);
end

function Y = kron_it(T,X)
% this function computes kron(I,A)*x
Y=T*X;
Y=Y(:);
% m =size(T,1);
% test=zeros(m^2,1);
% x=X(:);
% for i=1:m
%     test((i-1)*m+1 : i*m) = T*x((i-1)*m+1:i*m);
% end
% max(max(abs(Y-test)))
end

function Y = kron_ti(T,X)
% this function computes kron(A,I)*x
m=size(T,1);
m2=m^2;
Y=zeros(m2,1);
% for i=1:m
%     for j=1:m
%         Y((i-1)*m+j) = T(i,:)*X(j:m:m2);
%     end
% end
tmp=(0:m-1)*m;
for j=1:m
    Y(tmp+j) = T*X(j:m:m2);
end
end

