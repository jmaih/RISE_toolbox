function [T1,T0_T1]=functional_iteration_h_full(T0,Gplus01,A0,Aminus)
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


[~,n,h]=size(A0);

if size(Gplus01,2)~=n
    
    error('number of columns of Gplus01 inconsistent with the number of variables')
    
end

if size(Aminus,2)~=n
    
    error('number of columns of Aminus inconsistent with the number of variables')

end

T0=reshape(T0,[n,n,h]);

T1=T0;

for r0=1:h
    
    U=A0(:,:,r0);
    
    for r1=1:h
        
        U=U+Gplus01(:,:,r0,r1)*T0(:,:,r1);
    
    end
    
    T1(:,:,r0)=-U\Aminus(:,:,r0);

end

% update T
%---------
T1=reshape(T1,[n,n*h]);

T0_T1=reshape(T0,[n,n*h])-T1;

end