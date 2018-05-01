function [H,iH,dH] = conditioner(H0,tol)
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

% based on ddom.m
% returns a diagonally dominant matrix H by modifying the
% diagonal of H0.

if nargin<2
    tol=100*eps;
end
[m,n] = size(H0);
if m~=n
    error([mfilename,':: matrix must be square'])
end

if any(any(isnan(H0)))||any(any(isinf(H0)))
    H=H0;
    iH=nan(size(H));
    dH=nan;
else
    H0=make_diagonally_dominant(H0);
    [v,d] = eig(H0);
    d=diag(d);
    bad=d<0;
    if any(bad)
        d(bad)=tol;
    end
    iH=v*diag(1./d)*v';
    H=v*diag(d)*v';
    dH=prod(d);
end


% H=H0;
% if loop
%     for ii=1:n
%         d=H0(ii,ii);
%         a=abs(d);
%         f=0;
%         for jj=1:n
%             if ii~=jj
%                 f=f+abs(H0(ii,jj));
%             end
%         end
%         if f>=a
%             aii=(1+tol)*max(f,tol);
%             if d<0
%                 aii=-aii;
%             end
%             H(ii,ii)=aii;
%         end
%     end
% else
%     d = diag(H0);
%     a = abs(d);
%     f = sum(abs(H0),2)-a;
%     ii = find(f >= a);
%     k = ii + (ii-1)*m;
%     s = 2*(d(ii)>=0)-1;
%     H(k) = (1+tol)*s.*max(f(ii),tol);
% end

    function H=make_diagonally_dominant(H0)
        H=H0;
        d = diag(H0);
        a = abs(d);
        f = sum(abs(H0),2)-a;
        ii = find(f >= a);
        k = ii + (ii-1)*m;
        s = 2*(d(ii)>=0)-1;
        H(k) = (1+tol)*s.*max(f(ii),tol);
    end
end
