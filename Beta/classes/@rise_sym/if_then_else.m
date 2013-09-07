function c=if_then_else(state,a,b)
if nargin<3
	b=0;
end
	
if ~isa(a,'rise_sym')
    a=rise_sym(a);
elseif ~isa(b,'rise_sym')
    b=rise_sym(b); 
end
   c=rise_sym.multinary_operation('if_then_else',state,a,b);
end
