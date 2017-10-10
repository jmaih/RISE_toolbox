function [T,R,eigvals]=companion(B,ShockImpact,nz)

% the VBR has the form y=C*z+B1*y{-1}+B2*y{-2}+...+Bp*y{-p}+ShockImpact*u
% so that B=[C,B1,B2,...,Bp];
% nz: number of deterministic terms

B=B(:,nz+1:end);

[n,np]=size(B);

p=np/n;

dd=n*(p-1);

T=[B;
    eye(dd),zeros(dd,n)];

R=[];

if ~isempty(ShockImpact)

	R=[ShockImpact;zeros(dd,n)];
	
end

if nargout>2
    
    eigvals=abs(eig(T));
    
end

end
