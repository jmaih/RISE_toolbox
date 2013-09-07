function c=if_then_else(state,a,b)
if nargin<3
	b=0;
end
	
if ~isa(a,'planar')
    a=planar(a);
elseif ~isa(b,'planar')
    b=planar(b); 
end
   c=planar.multinary_operation('if_then_else',state,a,b);
end
