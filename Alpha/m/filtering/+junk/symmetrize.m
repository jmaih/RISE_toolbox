function B=symmetrize(A)%,flag
% if nargin<2
%     flag=1;
% end
% switch flag
%     case 1
        B=.5*(A+A');
%     otherwise
%         B=triu(A);
%         B=B+triu(B,1)';
% end
