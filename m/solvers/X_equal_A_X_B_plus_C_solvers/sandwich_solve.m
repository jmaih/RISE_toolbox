function [X,retcode]=sandwich_solve(A,B,C)
% sandwich_solve solves the linear equation X=A*X*B+C
%
% ::
%
%   [X,retcode]=sandwich_solve(A,B,C)
%
% Args:
%    - A :
%    - B :
%    - C :
%
% Returns:
%    :
%    - X :
%    - retcode :
%
% Note:
%
% Example:
%
%    See also:

% 
[U,R]=schur(B,'complex'); 
% we take the complex just to make sure that R is upper triangular
n=size(A,1);
I=eye(n);
Y=nan(n);
AY=nan(n);
F=C*U;
retcode=0;
for k=1:n
    tmp=(I-R(k,k)*A);
    if rcond(tmp)<1e-12
        X=nan(n);
        retcode=302;
        return
    end
    Y(:,k)=tmp\(AY(:,1:k-1)*R(1:k-1,k)+F(:,k));
    if k<n
        % in order to avoid recomputing AY all the time.
        AY(:,k)=A*Y(:,k);
    end
end
% although the imaginary part is potentially very tiny, we have to real
% everything since we took a complex decomposition above.
X=real(Y*U');

% % function [E,retcode]=sandwich_solve(A,B,C,hessenberg)
% % % solves the problem E=A*E*B+C
% % % for low dimensions, the result can be computed as
% % % E=reshape((eye(n^2)-kron(B',A))\C(:),n,n)
% % 
% % if nargin<4
% %     hessenberg=false;
% %     if nargin<3
% %         error([mfilename,':: need 3 or 4 arguments'])
% %     end
% % elseif nargin>4
% %     error([mfilename,':: need 3 or 4 arguments'])
% % end
% % 
% % retcode=0;
% % if any(any(isnan(A)))||any(any(isnan(B)))||any(any(isnan(C)))||...
% %         any(any(isinf(A)))||any(any(isinf(B)))||any(any(isinf(C)))
% %     retcode=1;
% %     E=[];
% %     return
% % end
% % n=size(A,1);
% % % schur decomposition of terms premultiplying E
% % I=eye(n);
% % if hessenberg
% %     [T,S,W,Z] = hess(I,A);
% % else
% %     [T,S,W,Z] = qz(I,A,'real');
% % end
% % % schur decomposition of the term post-multiplying E
% % [U,R] =schur(B,'complex');
% % % here we need to take the complex schur decomposition to ensure that R is
% % % upper triangular. The real schur of matlab only delivers a
% % % quasi-triangualar matrix.
% % 
% % F=W*C*U;
% % Y=nan(n);
% % for k=1:n
% %     Y(:,k)=(T-R(k,k)*S)\(S*Y(:,1:k-1)*R(1:k-1,k)+F(:,k)) ;
% % end
% % % SY=nan(n);
% % % for k=1:n
% % %     Y(:,k)=(T-R(k,k)*S)\(SY(:,1:k-1)*R(1:k-1,k)+F(:,k));
% % %     SY(:,k)=S*Y(:,k);
% % % end
% % 
% % % we've taken the complex decomposition of R above. Now we amend for that
% % % by taking the real of the result.
% % E=real(Z*Y*U');
% % 
% % % test=max(max(abs(E-(A*E*B+C))));
% % if any(any(isnan(E)))||any(any(isinf(E)))% ||test>1e-6
% %     retcode=1;
% % end;
% % 
% % % the algebra is as follows
% % % E=A*E*B+C ==> W*E*U=W*(A*E*B)*U+W*C*U
% % % ==> W*Z*Z'*E*U=W*A*Z*Z'*E*B*U+W*C*U
% % % ==> T*Z'*E*U=S*Z'*E*U*U'*B*U+F
% % % ==> T*Y=S*Y*U'*B*U+F
% % % ==> T*Y=S*Y*R+F
% % % ==> T*Y=S*[Y(:,1),...,Y(:,k),...,Y(:,n)]*R+F
% % % equating the kth column
% % % ==> T*Y(:,k)=S*[Y(:,1),...,Y(:,k),...,Y(:,n)]*R(:,k)+F(:,k)
% % % ==> T*Y(:,k)=S*[Y(:,1)*R(1,k),...,Y(:,k)*R(k,k),...,Y(:,n)*R(n,k)]+F(:,k)
% % % exploiting the fact that R is upper triangular
% % % ==> T*Y(:,k)-S*Y(:,k)*R(k,k)=S*[Y(:,1)*R(1,k),...,Y(:,k-1)*R(k-1,k)]+F(:,k)
% % % ==> [T-R(k,k)*S]*Y(:,k)=S*[Y(:,1)*R(1,k),...,Y(:,k-1)*R(k-1,k)]+F(:,k)
% % % ==> Y(:,k)=[T-R(k,k)*S]\(S*[Y(:,1)*R(1,k),...,Y(:,k-1)*R(k-1,k)]+F(:,k))
% % % ==> Y(:,k)=[T-R(k,k)*S]\(S*Y(:,1:k-1)*R(1:k-1,k)+F(:,k))
