function ypred=locate_permutation(x,n,wrt_wise)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

if nargin<3
    % [111 211 311 121 221 321 131 231 331 112 212 312 122 222 322 132 232
    % 332 113 213 313 123 223 323 133 233 333] <== wrt_wise
    wrt_wise=false;
    % [111 112 113 121 122 123 131 132 133 211 212 213 221 222 223 231 232
    % 233 311 312 313 321 322 323 331 332 333] <== ~wrt_wise
end
xlong=x(:);
if any(xlong>n)||any(xlong<=0)||any(floor(xlong)~=ceil(xlong))
    error('wrong specification of the elements in the x matrix')
end
order=size(x,2);
if wrt_wise
    ypred=x(:,end)-1;
    for ii=order-1:-1:1
        ypred=n*ypred+x(:,ii)-1;
    end
    ypred=ypred+1;
else
    ypred=x(:,1)-1;
    for ii=2:order
        ypred=n*ypred+x(:,ii)-1;
    end
    ypred=ypred+1;
end
end