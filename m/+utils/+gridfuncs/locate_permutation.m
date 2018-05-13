function ypred=locate_permutation(x,n,wrt_wise)
% locate_permutation locates a permutation deriving from differentiation
%
% ::
%
%
%   ypred=locate_permutation(x,n,wrt_wise)
%
% Args:
%
%    - **x** [scalar,vector|matrix]: nperm x order matrix of the permutations
%      of interest. The number of rows represents the number of permutations
%      and the number of columns the order of differentiation
%
%    - **n** [scalar]: number of variables in the differentiation
%
%    - **wrt_wise** [true|{false}]: when false, the retrieved order is that of
%      a kronecker unfolded column wise i.e. [111 112 113 121 122 123 131 132
%      133 211 212 213 221 222 223 231 232 233 311 312 313 321 322 323 331 332
%      333]. Whe true, the order is [111 211 311 121 221 321 131 231 331 112
%      212 312 122 222 322 132 232 332 113 213 313 123 223 323 133 233 333]
%
% Returns:
%    :
%
%    - **ypred** [scalar|vector] : location of the permutations
%
% Note:
%
% Example:
%
%    See also:

if nargin<3

    wrt_wise=false;

end

% xlong=x(:); % checking this slows things down
% 
% if any(xlong>n)||any(xlong<=0)||any(floor(xlong)~=ceil(xlong))
% 
%     error('wrong specification of the elements in the x matrix')
% 
% end

order=size(x,2);

if wrt_wise
    
    ypred=x(:,end)-1;
    
    for ii=order-1:-1:1
    
        ypred=n*ypred+x(:,ii)-1;
    
    end
    
else
    
    ypred=x(:,1)-1;
    
    for ii=2:order
    
        ypred=n*ypred+x(:,ii)-1;
    
    end

end
    
ypred=ypred+1;

end