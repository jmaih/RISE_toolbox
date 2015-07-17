function c=kron(a,b)
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
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

a=splanar(a);
b=splanar(b);

c=splanar.empty(0);
[ra,ca]=size(a);
[rb,cb]=size(b);

irow=0;
for ia=1:ra
	for ib=1:rb
		irow=irow+1;
		jcol=0;
		for ja=1:ca
			for jb=1:cb
				jcol=jcol+1;
				c(irow,jcol)=a(ia,ja)*b(ib,jb);
			end
		end
	end
end
