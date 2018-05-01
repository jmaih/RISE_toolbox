classdef aplanar
% aplanar automatic "planar" differentiation
%
% - [abs](aplanar/abs)
% - [acos](aplanar/acos)
% - [acosh](aplanar/acosh)
% - [and](aplanar/and)
% - [aplanar](aplanar/aplanar)
% - [asin](aplanar/asin)
% - [asinh](aplanar/asinh)
% - [atan](aplanar/atan)
% - [atanh](aplanar/atanh)
% - [cos](aplanar/cos)
% - [cosh](aplanar/cosh)
% - [cot](aplanar/cot)
% - [coth](aplanar/coth)
% - [diff](aplanar/diff)
% - [eq](aplanar/eq)
% - [erf](aplanar/erf)
% - [exp](aplanar/exp)
% - [ge](aplanar/ge)
% - [gt](aplanar/gt)
% - [if_elseif](aplanar/if_elseif)
% - [if_then_else](aplanar/if_then_else)
% - [le](aplanar/le)
% - [log](aplanar/log)
% - [log10](aplanar/log10)
% - [lt](aplanar/lt)
% - [max](aplanar/max)
% - [min](aplanar/min)
% - [minus](aplanar/minus)
% - [mpower](aplanar/mpower)
% - [mrdivide](aplanar/mrdivide)
% - [mtimes](aplanar/mtimes)
% - [ne](aplanar/ne)
% - [normcdf](aplanar/normcdf)
% - [normpdf](aplanar/normpdf)
% - [or](aplanar/or)
% - [plus](aplanar/plus)
% - [power](aplanar/power)
% - [rdivide](aplanar/rdivide)
% - [sign](aplanar/sign)
% - [sin](aplanar/sin)
% - [sinh](aplanar/sinh)
% - [sqrt](aplanar/sqrt)
% - [tan](aplanar/tan)
% - [tanh](aplanar/tanh)
% - [times](aplanar/times)
% - [uminus](aplanar/uminus)
% - [uplus](aplanar/uplus)
properties
x
dx
dxx
dxxx
dxxxx
dxxxxx
order
nvars
deriv_names
end
methods
% constructor
%------------
function obj=aplanar(v,nvars,vpos,order_)
% aplanar --  constructor for aplanar differentiation
%
% ::
%
%
%   obj=aplanar(v,nvars,vpos,order_)
%
% Args:
%              %
%              % - **v** [vector]:
%              %
%              % - **nvars** [integer]:
%              %
%              % - **vpos** [scalar|vector]:
%              %
%              % - **order_** [integer]:
%              %
% Returns:
%    :
%              %
%              % - **obj** [aplanar]: aplanar object with derivatives
%              %
% Note:
%              %
%              % - it is desirable to further save some computations by
%              % computing the shrinks only and then expanding entire maps as
%              % in splanar...
%              %
% Example:
%
% See also:
if nargin
nv=numel(v);
if nv~=numel(vpos)
error('number of positions should be the same as the number of locations')
end
if order_<=0||(order_~=floor(order_))
error('order of differentiation should be a strictly positive integer')
end
obj.order=order_;
obj.nvars=nvars;
further=repmat({nvars},1,order_);
xxx=repmat('x',1,order_);
dnames=cell(1,order_);
for io=1:order_
dd=['d',xxx(1:io)];
ff=further(1:io);
obj.(dd)=zeros(1,ff{:});
dnames{io}=dd;
end
obj.deriv_names=dnames;
if ~isnan(vpos(1))
obj.dx(vpos(1))=1;
end
obj.x=v(1);
if nv>1
tmp=obj;
for ii=2:nv
obj.x=v(ii);
obj.dx=zeros(1,nvars);
if ~isnan(vpos(ii))
obj.dx(vpos(ii))=1;
end
tmp(ii,1)=obj;
end
obj=tmp;
end
end
end
% unary
%------
function obj = abs(g)
obj=apply_to_all(g,@times,sign(g.x));
end

function obj = acos(g)
obj=0.5*pi-asin(g);
end

function obj = acosh(g)
obj=2*log(sqrt(.5*(g+1))+sqrt(.5*(g-1)));
end

function obj = asin(g)
obj=-1i*log(1i*g+sqrt(1-g^2));
end

function obj = asinh(g)
obj=log(g+sqrt(1+g^2));
end

function obj = atan(g)
obj=-1i*log((1+1i*g)*sqrt(1/(1+g^2)));
end

function obj = atanh(g)
obj=.5*(log(1+g)-log(1-g));
end

function obj = cos(g)
obj=sine_cosine(g,'cos');
end

function obj = cosh(g)
obj=.5*(exp(g)+exp(-g));
end

function obj = cot(g)
obj=1/tan(g);
end

function obj = coth(g)
obj = cosh(g)/sinh(g);
end

function obj = erf(g)
obj=reinitialize(g);
f=erf(g.x);
obj.x=f;
oo=obj.order;
for ii=1:g.nvars
if ii==1
obj.dx=2/sqrt(pi)*exp(-g.x^2)*g.dx;
if oo==1
break
end
end
if oo>1 && obj.dx(1,ii)~=0
for jj=1:ii
do_second_order();
if oo>2 && obj.dxx(1,ii,jj)~=0
for kk=1:jj
do_third_order();
if oo>3 && obj.dxxx(1,ii,jj,kk)~=0
for ll=1:kk
do_fourth_order();
if oo>4 && obj.dxxxx(1,ii,jj,kk,ll)~=0
for mm=1:ll
do_fifth_order();
end
end
end
end
end
end
end
end
end

function do_second_order()
obj.dxx(1,ii,jj)=2*(g.dxx(1,ii,jj)*f-g.dx(1,jj)*obj.dx(1,ii));
end

function do_third_order()
obj.dxxx(1,ii,jj,kk)=2*(g.dxxx(1,ii,jj,kk)*f...
+g.dxx(1,ii,jj)*obj.dx(1,kk)...
-g.dxx(1,jj,kk)*obj.dx(1,ii)...
-g.dx(1,jj)*g.dxx(1,ii,kk));
end

function do_fourth_order()
obj.dxxxx(1,ii,jj,kk,ll)=2*(...
g.dxxx(1,ii,jj,kk,ll)*f...
+g.dxxx(1,ii,jj,kk)*obj.dx(1,ll)...
+g.dxxx(1,ii,jj,ll)*obj.dx(1,kk)...
+g.dxx(1,ii,jj)*obj.dxx(1,kk,ll)...
-g.dxxx(1,jj,kk,ll)*obj.dx(1,ii)...
-g.dxx(1,jj,kk)*obj.dxx(1,ii,ll)...
-g.dxx(1,jj,ll)*obj.dxx(1,ii,kk)...
-g.dx(1,jj)*obj.dxxx(1,ii,kk,ll)...
);
end

function do_fifth_order()
obj.dxxxxx(1,ii,jj,kk,ll,mm)=2*(...
g.dxxxx(1,ii,jj,kk,ll,mm)*f...
+g.dxxxx(1,ii,jj,kk,ll)*obj.dx(1,mm)...
+g.dxxxx(1,ii,jj,kk,mm)*obj.dx(1,ll)...
+g.dxxx(1,ii,jj,kk)*obj.dxx(1,ll,mm)...
+g.dxxxx(1,ii,jj,ll,mm)*obj.dx(1,kk)...
+g.dxxx(1,ii,jj,ll)*obj.dxx(1,kk,mm)...
+g.dxxx(1,ii,jj,mm)*obj.dxx(1,kk,ll)...
+g.dxx(1,ii,jj)*obj.dxxx(1,kk,ll,mm)...
-g.dxxxx(1,jj,kk,ll,mm)*obj.dx(1,ii)...
-g.dxxx(1,jj,kk,ll)*obj.dxx(1,ii,mm)...
-g.dxxx(1,jj,kk,mm)*obj.dxx(1,ii,ll)...
-g.dxx(1,jj,kk)*obj.dxxx(1,ii,ll,mm)...
-g.dxxx(1,jj,ll,mm)*obj.dxx(1,ii,kk)...
-g.dxx(1,jj,ll)*obj.dxxx(1,ii,kk,mm)...
-g.dxx(1,jj,mm)*obj.dxxx(1,ii,kk,ll)...
-g.dx(1,jj)*obj.dxxxx(1,ii,kk,ll,mm)...
);
end

end

function obj = exp(g)
obj=reinitialize(g);
f=exp(g.x);
obj.x=f;
oo=obj.order;
for ii=1:g.nvars
if ii==1
obj.dx=g.dx*f;
if oo==1
break
end
end
if oo>1 && obj.dx(1,ii)~=0
for jj=1:ii
do_second_order();
if oo>2 && obj.dxx(1,ii,jj)~=0
for kk=1:jj
do_third_order();
if oo>3 && obj.dxxx(1,ii,jj,kk)~=0
for ll=1:kk
do_fourth_order();
if oo>4 && obj.dxxxx(1,ii,jj,kk,ll)~=0
for mm=1:ll
do_fifth_order();
end
end
end
end
end
end
end
end
end

function do_second_order()
obj.dxx(1,ii,jj)=g.dxx(1,ii,jj)*f+g.dx(1,ii)*obj.dx(1,jj);
end

function do_third_order()
obj.dxxx(1,ii,jj,kk)=g.dxxx(1,ii,jj,kk)*f+...
g.dxx(1,ii,jj)*obj.dx(1,kk)+...
g.dxx(1,ii,kk)*obj.dx(1,jj)+...
g.dx(1,ii)*obj.dxx(1,jj,kk);
end

function do_fourth_order()
obj.dxxxx(1,ii,jj,kk,ll)=g.dxxxx(1,ii,jj,kk,ll)*f+...
g.dxxx(1,ii,jj,kk)*obj.dx(1,ll)+...
g.dxxx(1,ii,kk,ll)*obj.dx(1,jj)+...
g.dxx(1,ii,kk)*obj.dxx(1,jj,ll)+...
g.dxxx(1,jj,kk,ll)*obj.dx(1,ii)+...
g.dxx(1,jj,kk)*obj.dxx(1,ii,ll)+...
g.dxx(1,kk,ll)*obj.dxx(1,ii,jj)+...
g.dx(1,kk)*obj.dxxx(1,ii,jj,ll);
end

function do_fifth_order()
obj.dxxxxx(1,ii,jj,kk,ll,mm)=g.dxxxxx(1,ii,jj,kk,ll,mm)*f+...
g.dxxxx(1,ii,jj,kk,ll)*obj.dx(1,mm)+...
g.dxxxx(1,ii,jj,kk,mm)*obj.dx(1,ll)+...
g.dxxxx(1,ii,kk,ll,mm)*obj.dx(1,jj)+...
g.dxxxx(1,jj,kk,ll,mm)*obj.dx(1,ii)+...
g.dxxx(1,ii,jj,kk)*obj.dxx(1,ll,mm)+...
g.dxxx(1,ii,kk,ll)*obj.dxx(1,jj,mm)+...
g.dxxx(1,ii,kk,mm)*obj.dxx(1,jj,ll)+...
g.dxxx(1,jj,kk,ll)*obj.dxx(1,ii,mm)+...
g.dxxx(1,jj,kk,mm)*obj.dxx(1,ii,ll)+...
g.dxxx(1,kk,ll,mm)*obj.dxx(1,ii,jj)+...
g.dxx(1,ii,kk)*obj.dxxx(1,jj,ll,mm)+...
g.dxx(1,jj,kk)*obj.dxxx(1,ii,ll,mm)+...
g.dxx(1,kk,ll)*obj.dxxx(1,ii,jj,mm)+...
g.dxx(1,kk,mm)*obj.dxxx(1,ii,jj,ll)+...
g.dx(1,kk)*obj.dxxxx(1,ii,jj,ll,mm);
end

end

function obj = log(g)
obj=reinitialize(g);
obj.x=log(g.x);
oo=obj.order;
igx=1/g.x;
for ii=1:g.nvars
if ii==1
obj.dx=g.dx*igx;
if oo==1
break
end
end
if oo>1 && obj.dx(1,ii)~=0
for jj=1:ii
do_second_order();
if oo>2 && obj.dxx(1,ii,jj)~=0
for kk=1:jj
do_third_order();
if oo>3 && obj.dxxx(1,ii,jj,kk)~=0
for ll=1:kk
do_fourth_order();
if oo>4 && obj.dxxxx(1,ii,jj,kk,ll)~=0
for mm=1:ll
do_fifth_order();
end
end
end
end
end
end
end
end
end

function do_second_order()
obj.dxx(1,ii,jj)=igx*g.dxx(1,ii,jj)-obj.dx(1,ii)*obj.dx(1,jj);
end

function do_third_order()
obj.dxxx(1,ii,jj,kk)=igx*(g.dxxx(1,ii,jj,kk)-g.dxx(1,ii,jj)*obj.dx(1,kk))...
-(...
obj.dxx(1,jj,kk)*obj.dx(1,ii)...
+obj.dxx(1,ii,kk)*obj.dx(1,jj)...
);
end

function do_fourth_order()
obj.dxxxx(1,ii,jj,kk,ll)=igx*(...
g.dxxxx(1,ii,jj,kk,ll)...
-g.dxxx(1,ii,jj,kk)*obj.dx(1,ll)...
-(g.dxxx(1,ii,jj,ll)-g.dxx(1,ii,jj)*obj.dx(1,ll))*obj.dx(1,kk)...
-g.dxx(1,ii,jj)*obj.dxx(1,kk,ll)...
)-(...
obj.dxxx(1,jj,kk,ll)*obj.dx(1,ii)...
+obj.dxxx(1,ii,kk,ll)*obj.dx(1,jj)...
+obj.dxx(1,jj,kk)*obj.dxx(1,ii,ll)...
+obj.dxx(1,jj,ll)*obj.dxx(1,ii,kk));
end

function do_fifth_order()
a1=g.dxxxxx(1,ii,jj,kk,ll,mm)-g.dxxxx(1,ii,jj,kk,ll)*obj.dx(1,mm);

a2=(g.dxxxx(1,ii,jj,kk,mm)-g.dxxx(1,ii,jj,kk)*obj.dx(1,mm))*obj.dx(1,ll)+...
g.dxxx(1,ii,jj,kk)*obj.dxx(1,ll,mm);

a3=(g.dxxxx(1,ii,jj,ll,mm)-g.dxxx(1,ii,jj,ll)*obj.dx(1,mm))*obj.dx(1,kk)+...
g.dxxx(1,ii,jj,ll)*obj.dxx(1,kk,mm);

a4=(g.dxxx(1,ii,jj,mm)-g.dxx(1,ii,jj)*obj.dx(1,mm))*obj.dx(1,ll)*obj.dx(1,kk)+...
g.dxx(1,ii,jj)*obj.dxx(1,ll,mm)*obj.dx(1,kk)+...
g.dxx(1,ii,jj)*obj.dx(1,ll)*obj.dxx(1,kk,mm);

a5=(g.dxxx(1,ii,jj,mm)-g.dxx(1,ii,jj)*obj.dx(1,mm))*obj.dxx(1,kk,ll)+...
g.dxx(1,ii,jj)*obj.dxxx(1,kk,ll,mm);

a6=obj.dxxxx(1,ii,kk,ll,mm)*obj.dx(1,jj)+...
obj.dxxx(1,ii,kk,ll)*obj.dxx(1,jj,mm)+...
obj.dxxx(1,ii,kk,mm)*obj.dxx(1,jj,ll)+...
obj.dxx(1,ii,kk)*obj.dxxx(1,jj,ll,mm)+...
obj.dxxx(1,ii,ll,mm)*obj.dxx(1,jj,kk)+...
obj.dxx(1,ii,ll)*obj.dxxx(1,jj,kk,mm)+...
obj.dxx(1,ii,mm)*obj.dxxx(1,jj,kk,ll)+...
obj.dx(1,ii)*obj.dxxxx(1,jj,kk,ll,mm);

obj.dxxxxx(1,ii,jj,kk,ll,mm)=igx*(a1-a2-a3+a4-a5)-a6;
end

end

function obj = log10(g)
obj=log(g)/log(10);
end

function obj = sign(g)
obj=zero_derivatives(@sign,g);
end

function obj = sin(g)
obj=sine_cosine(g,'sin');
end

function obj = sinh(g)
obj=.5*(exp(g)-exp(-g));
end

function obj = sqrt(g)
obj=g^0.5;
end

function obj = steady_state(g)
obj=zero_derivatives(@steady_state,g);
end

function obj = tan(g)
obj = sin(g)/cos(g);
end

function obj = tanh(g)
obj = sinh(g)/cosh(g);
end

function obj = uminus(g)
obj=g;
index=obj.deriv_names;
obj.x=-g.x;
for ii=1:obj.order
ff=index{ii};
obj.(ff)=-g.(ff);
end
end

function obj = uplus(g)
obj=g;
end

% binary
%-------
function obj = and(g,h)
obj=zero_derivatives(@and,g,h);
end

function obj = eq(g,h)
obj=zero_derivatives(@eq,g,h);
end

function obj = ge(g,h)
obj=zero_derivatives(@ge,g,h);
end

function obj = gt(g,h)
obj=zero_derivatives(@gt,g,h);
end

function obj = le(g,h)
obj=zero_derivatives(@le,g,h);
end

function obj = lt(g,h)
obj=zero_derivatives(@lt,g,h);
end

function obj = max(g,h)
g_dbl=isa(g,'double');
h_dbl=isa(h,'double');
if g_dbl
obj=h;
if g>=obj.x
obj.x=g;
obj=zero_derivatives(@(x)x,obj);
end
elseif h_dbl
obj=g;
if h>=obj.x
obj.x=h;
obj=zero_derivatives(@(x)x,obj);
end
else
if h.x>=g.x
obj=h;
else
obj=g;
end
end
end

function obj = min(g,h)
obj = -max(-g,-h);
end

function obj = minus(g,h)
g_dbl=isa(g,'double');
h_dbl=isa(h,'double');
if g_dbl
obj=-h;
obj.x=obj.x+g;
elseif h_dbl
obj=g;
obj.x=g.x-h;
else
obj=g;
index=obj.deriv_names;
obj.x=g.x-h.x;
for ii=1:obj.order
ff=index{ii};
obj.(ff)=g.(ff)-h.(ff);
end
end
end

function obj = mpower(g,h)
small=1e-32;
g_dbl=isa(g,'double');
if ~g_dbl
if g.x==0
g.x=small;
end
elseif g_dbl && g==0
g=small;
end
obj=exp(h*log(g));
end

function obj = mrdivide(g,h)
obj=g*h^(-1);
end

function obj = mtimes(g,h)
g_dbl=isa(g,'double');
h_dbl=isa(h,'double');
if g_dbl
% if g==0,obj=0;
obj=apply_to_all(h,@times,g);
elseif h_dbl
% if h==0,obj=0;
obj=apply_to_all(g,@times,h);
else
obj=reinitialize(g);
obj.x=g.x*h.x;
oo=obj.order;
for ii=1:g.nvars
if ii==1
obj.dx=g.x*h.dx+h.x*g.dx;
if oo==1
break
end
end
if oo>1 && obj.dx(1,ii)~=0
for jj=1:ii
do_second_order();
if oo>2 && obj.dxx(1,ii,jj)~=0
for kk=1:jj
do_third_order();
if oo>3 && obj.dxxx(1,ii,jj,kk)~=0
for ll=1:kk
do_fourth_order();
if oo>4 && obj.dxxxx(1,ii,jj,kk,ll)~=0
for mm=1:ll
do_fifth_order();
end
end
end
end
end
end
end
end
end
end

function do_second_order()
obj.dxx(1,ii,jj)=g.dxx(1,ii,jj)*h.x...
+g.dx(1,ii)*h.dx(1,jj)...
+g.dx(1,jj)*h.dx(1,ii)...
+g.x*h.dxx(1,ii,jj);
end

function do_third_order()
obj.dxxx(1,ii,jj,kk)=g.dxxx(1,ii,jj,kk)*h.x...
+g.dxx(1,ii,jj)*h.dx(1,kk)...
+g.dxx(1,ii,kk)*h.dx(1,jj)...
+g.dxx(1,jj,kk)*h.dx(1,ii)...
+g.dx(1,ii)*h.dxx(1,jj,kk)...
+g.dx(1,jj)*h.dxx(1,ii,kk)...
+g.dx(1,kk)*h.dxx(1,ii,jj)...
+g.x*h.dxxx(1,ii,jj,kk);
end

function do_fourth_order()
obj.dxxxx(1,ii,jj,kk,ll)=g.dxxxx(1,ii,jj,kk,ll)*h.x...
+g.dxxx(1,ii,jj,kk)*h.dx(1,ll)...
+g.dxxx(1,ii,jj,ll)*h.dx(1,kk)...
+g.dxxx(1,ii,kk,ll)*h.dx(1,jj)...
+g.dxxx(1,jj,kk,ll)*h.dx(1,ii)...
+g.dxx(1,ii,jj)*h.dxx(1,kk,ll)...
+g.dxx(1,ii,kk)*h.dxx(1,jj,ll)...
+g.dxx(1,ii,ll)*h.dxx(1,jj,kk)...
+g.dxx(1,jj,kk)*h.dxx(1,ii,ll)...
+g.dxx(1,jj,ll)*h.dxx(1,ii,kk)...
+g.dxx(1,kk,ll)*h.dxx(1,ii,jj)...
+g.dx(1,jj)*h.dxxx(1,ii,kk,ll)...
+g.dx(1,kk)*h.dxxx(1,ii,jj,ll)...
+g.dx(1,ll)*h.dxxx(1,ii,jj,kk)...
+g.dx(1,ii)*h.dxxx(1,jj,kk,ll)...
+g.x*h.dxxxx(1,ii,jj,kk,ll);
end

function do_fifth_order()
obj.dxxxxx(1,ii,jj,kk,ll,mm)=g.dxxxxx(1,ii,jj,kk,ll,mm)*h.x...
+g.dxxxx(1,ii,jj,kk,ll)*h.dx(1,mm)...
+g.dxxxx(1,ii,jj,kk,mm)*h.dx(1,ll)...
+g.dxxxx(1,ii,kk,ll,mm)*h.dx(1,jj)...
+g.dxxxx(1,jj,kk,ll,mm)*h.dx(1,ii)...
+g.dxxxx(1,ii,jj,ll,mm)*h.dx(1,kk)...
+g.dxxx(1,ii,jj,kk)*h.dxx(1,ll,mm)...
+g.dxxx(1,ii,jj,mm)*h.dxx(1,kk,ll)...
+g.dxxx(1,ii,kk,ll)*h.dxx(1,jj,mm)...
+g.dxxx(1,ii,kk,mm)*h.dxx(1,jj,ll)...
+g.dxxx(1,ii,jj,ll)*h.dxx(1,kk,mm)...
+g.dxxx(1,ii,ll,mm)*h.dxx(1,jj,kk)...
+g.dxxx(1,jj,kk,mm)*h.dxx(1,ii,ll)...
+g.dxxx(1,jj,ll,mm)*h.dxx(1,ii,kk)...
+g.dxxx(1,kk,ll,mm)*h.dxx(1,ii,jj)...
+g.dxxx(1,jj,kk,ll)*h.dxx(1,ii,mm)...
+g.dxx(1,ii,jj)*h.dxxx(1,kk,ll,mm)...
+g.dxx(1,ii,kk)*h.dxxx(1,jj,ll,mm)...
+g.dxx(1,ii,ll)*h.dxxx(1,jj,kk,mm)...
+g.dxx(1,ii,mm)*h.dxxx(1,jj,kk,ll)...
+g.dxx(1,jj,kk)*h.dxxx(1,ii,ll,mm)...
+g.dxx(1,jj,ll)*h.dxxx(1,ii,kk,mm)...
+g.dxx(1,jj,mm)*h.dxxx(1,ii,kk,ll)...
+g.dxx(1,kk,ll)*h.dxxx(1,ii,jj,mm)...
+g.dxx(1,kk,mm)*h.dxxx(1,ii,jj,ll)...
+g.dxx(1,ll,mm)*h.dxxx(1,ii,jj,kk)...
+g.dx(1,ii)*h.dxxxx(1,jj,kk,ll,mm)...
+g.dx(1,jj)*h.dxxxx(1,ii,kk,ll,mm)...
+g.dx(1,kk)*h.dxxxx(1,ii,jj,ll,mm)...
+g.dx(1,ll)*h.dxxxx(1,ii,jj,kk,mm)...
+g.dx(1,mm)*h.dxxxx(1,ii,jj,kk,ll)...
+g.x*h.dxxxxx(1,ii,jj,kk,ll,mm);
end

end

function obj = ne(g,h)
obj=zero_derivatives(@ne,g,h);
end

function obj = or(g,h)
obj=zero_derivatives(@or,g,h);
end

function obj = plus(g,h)
g_dbl=isa(g,'double');
h_dbl=isa(h,'double');
if g_dbl
obj=h;
obj.x=h.x+g;
elseif h_dbl
obj=g;
obj.x=g.x+h;
else
obj=g;
index=obj.deriv_names;
obj.x=g.x+h.x;
for ii=1:obj.order
ff=index{ii};
obj.(ff)=g.(ff)+h.(ff);
end
end
end

function obj = power(g,h)
obj=mpower(g,h);
end

function obj = rdivide(g,h)
obj=mrdivide(g,h);
end

function obj = times(g,h)
obj=mtimes(g,h);
end

% n-ary
%------
function obj = if_elseif(varargin)
% just pick the first element corresponding to the location
% that evaluates to true
nargs=length(varargin);
if rem(nargs,2)
error('the number of arguments must be even')
end
done=false;
iarg=1;
while ~done && iarg<nargs
check=varargin{iarg};
if ~(isa(check,'double')||isa(check,'logical'))
check=check.x;
end
if check
obj=varargin{iarg+1};
done=~done;
end
iarg=iarg+2;
end
if ~done
error('could not find a valid statement in if_elseif')
end
end

function obj = if_then_else(a,b,c)
obj = if_elseif(a,b,1-a,c);
end

function obj = normcdf(x,mu,sig)
if nargin<3
sig=[];
if nargin<2
mu=[];
end
end
if isempty(sig)
sig=1;
end
if isempty(mu)
mu=0;
end
obj=0.5*(1+erf((x-mu)/(sig*sqrt(2))));
end

function obj = normpdf(x,mu,sig)
if nargin<3
sig=[];
if nargin<2
mu=[];
end
end
if isempty(sig)
sig=1;
end
if isempty(mu)
mu=0;
end
obj=exp(-.5*((x-mu)/sig)^2)/(sig*sqrt(2*pi));
end

end

methods(Access=private)

function obj=apply_to_all(obj,func,const)
obj.x=func(obj.x,const);
index=obj.deriv_names;
for ii=1:obj.order
ff=index{ii};
obj.(ff)=func(obj.(ff),const);
end
end

function obj=zero_derivatives(func,a,b)
a_aplanar = isa(a,'aplanar');
if nargin==3
if a_aplanar
obj=a;
ax=a.x;
bx=b;
else
obj=b;
bx=b.x;
ax=a;
end
args={ax,bx};
elseif nargin==2
obj=a;
ax=a.x;
args={ax};
end
obj.x=func(args{:});
index=obj.deriv_names;
for ii=1:obj.order
obj.(index{ii})(:)=0;
end
end

end

methods(Static)

function C=diff(func,active,inactive,order)
n=numel(func);
proto_rows=cell(1,order);
nv=size(active,1);
recorder=struct('iter',{},'ncells',{},'info',{});
ncells=1000;
for io=1:order
proto_rows{io}=zeros([1,nv*ones(1,io)]);% C{ii}=zeros([n,nv*ones(1,ii)]);
recorder(io).iter=0;
recorder(io).ncells=ncells;
recorder(io).info=cell(ncells,1);
end
silent=true;
xxx=repmat('x',1,order);
proto_reorder=cell(1,order);
% bigset=cell(1,order);
for ifunc=1:n
% get info
%---------
objective=func{ifunc};
arg_list=get_arg_list();
% locate the variables for efficient differentiation
%---------------------------------------------------
pos_active=locate_variables(arg_list,active(:,1),silent);
if isempty(inactive)
pos_inactive=nan*pos_active;
else
pos_inactive=locate_variables(arg_list,inactive(:,1),silent);
end
nwrt=sum(~isnan(pos_active));
nargs=numel(arg_list);
vpos=0;
x0=[];
for iarg=1:nargs
if isnan(pos_active(iarg))
arg_list{iarg}=inactive{pos_inactive(iarg),2};
else
vpos=vpos+1;
if isempty(x0)
x0=aplanar(active{pos_active(iarg),2},nwrt,vpos,order);
else
x0.x=active{pos_active(iarg),2};
x0.dx(x0.dx~=0)=0;
x0.dx(vpos)=1;
end
arg_list{iarg}=x0;
end
end
% differentiate
%--------------
obj=objective(arg_list{:});
% store the information
%----------------------
pos_active(isnan(pos_active))=[];
if ~isempty(pos_active)
for io=1:order
further=repmat({pos_active},1,io);
tmp=proto_rows{io};
tmp(1,further{:})=redesign(obj.(['d',xxx(1:io)]));
tmp=reshape(tmp,1,[]);
nice=find(tmp);
record_these_derivatives(ifunc,nice,tmp(nice))
end
end
end
% now create output matrices
%----------------------------
max_bytes=utils.windows_mac.maxbytes()/8;
% MemAvailableAllArrays = MaxPossibleArrayBytes so we divide by
% 8 above to be economical...
C=cell(1,order);
for io=1:order
info=cell2mat(recorder(io).info(1:recorder(io).iter));
% minimum storage requirements: number of rows irrelevant
% 8 * nzmax + 8 * (nzmax + ncols + 1)
nzmax_=size(info,1);
ncols=nv^io;
mem_req=8*nzmax_+8*(nzmax_+ncols+1);
short_and_wide=mem_req<max_bytes;
tall_and_thin=~short_and_wide;
if short_and_wide
% store normally
C{io}=sparse(info(:,1),info(:,2),info(:,3),n,nv^io);
elseif tall_and_thin
% store and as transpose to save on memory
C{io}=tsparse(info(:,1),info(:,2),info(:,3),n,nv^io);
end
end

function record_these_derivatives(ifunc,pos,d)
iter=recorder(io).iter+1;
recorder(io).iter=recorder(io).iter+1;
if iter>=recorder(io).ncells
recorder(io).ncells=recorder(io).ncells+ncells;
recorder(io).info{end+ncells}={};
end
npos=numel(pos);
recorder(io).info{iter}=[ifunc(ones(npos,1)),pos(:),d(:)];
end

function C=redesign(C)
% assigns symmetric values to their location for orders
% greater than 1
if io>1
if isempty(proto_reorder{io})
proto_reorder{io}=cell2mat(utils.gridfuncs.mypermutation(1:io));
end
proto=proto_reorder{io};
C=autofill_symmetry(C,proto); % C=autofill_symmetry(C,proto,bigset(1:io),1,nwrt);
end
end
function arg_list=get_arg_list()
fstr=func2str(objective);
right_par=find(fstr==')',1,'first');
fstr=fstr(3:right_par-1);
arg_list= regexp(fstr,',','split');
end
end

end
end

function obj=sine_cosine(g,type)
switch type
case 'cos'
f=cos(g.x);
z=-sin(g.x);
case 'sin'
f=sin(g.x);
z=cos(g.x);
otherwise
error(['unknown type ',type])
end
obj=reinitialize(g);

obj.x=f;
oo=obj.order;
for ii=1:obj.nvars
if ii==1
obj.dx=g.dx*z;
if oo==1
break
end
end
if oo>1 && obj.dx(1,ii)~=0
for jj=1:ii
do_second_order();
if oo>2 && obj.dxx(1,ii,jj)~=0
for kk=1:jj
do_third_order();
if oo>3 && obj.dxxx(1,ii,jj,kk)~=0
for ll=1:kk
do_fourth_order();
if oo>4 && obj.dxxxx(1,ii,jj,kk,ll)~=0
for mm=1:ll
do_fifth_order();
end
end
end
end
end
end
end
end
end

function do_second_order()
obj.dxx(1,ii,jj)=z*g.dxx(1,ii,jj)-g.dx(1,ii)*g.dx(1,jj)*f;
end
function do_third_order()
obj.dxxx(1,ii,jj,kk)=-(g.dx(1,kk)*g.dxx(1,ii,jj)+...
g.dxx(1,jj,kk)*g.dx(1,ii)+...
g.dx(1,jj)*g.dxx(1,ii,kk))*f...
+z*g.dxxx(1,ii,jj,kk)...
-g.dx(1,ii)*g.dx(1,jj)*obj.dx(1,kk);
end
function do_fourth_order()
obj.dxxxx(1,ii,jj,kk,ll)=-(...
g.dxx(1,kk,ll)*g.dxx(1,ii,jj)+...
g.dx(1,kk)*g.dxxx(1,ii,jj,ll)+...
g.dxxx(1,jj,kk,ll)*g.dx(1,ii)+...
g.dxx(1,jj,kk)*g.dxx(1,ii,ll)+...
g.dxx(1,jj,ll)*g.dxx(1,ii,kk)+...
g.dx(1,jj)*g.dxxx(1,ii,kk,ll)+...
g.dx(1,ll)*g.dxxx(1,ii,jj,kk)...
)*f...
-(...
g.dx(1,kk)*g.dxx(1,ii,jj)+...
g.dxx(1,jj,kk)*g.dx(1,ii)+...
g.dx(1,jj)*g.dxx(1,ii,kk)...
)*obj.dx(1,ll)...
-(...
g.dxx(1,jj,ll)*g.dx(1,ii)+...
g.dx(1,jj)*g.dxx(1,ii,ll)...
)*obj.dx(1,kk)...
-g.dx(1,jj)*g.dx(1,ii)*obj.dxx(1,kk,ll)...
+z*g.dxxxx(1,ii,jj,kk,ll);
end
function do_fifth_order()
obj.dxxxxx(1,ii,jj,kk,ll,mm)=...
z*g.dxxxxx(1,ii,jj,kk,ll,mm)-(...
(...
g.dx(1,mm)*g.dxxxx(1,ii,jj,kk,ll)+...
g.dxxx(1,kk,ll,mm)*g.dxx(1,ii,jj)+g.dxx(1,kk,ll)*g.dxxx(1,ii,jj,mm)+...
g.dxx(1,kk,mm)*g.dxxx(1,ii,jj,ll)+g.dx(1,kk)*g.dxxxx(1,ii,jj,ll,mm)+...
g.dxxxx(1,jj,kk,ll,mm)*g.dx(1,ii)+g.dxxx(1,jj,kk,ll)*g.dxx(1,ii,mm)+...
g.dxxx(1,jj,kk,mm)*g.dxx(1,ii,ll)+g.dxx(1,jj,kk)*g.dxxx(1,ii,ll,mm)+...
g.dxxx(1,jj,ll,mm)*g.dxx(1,ii,kk)+g.dxx(1,jj,ll)*g.dxxx(1,ii,kk,mm)+...
g.dxx(1,jj,mm)*g.dxxx(1,ii,kk,ll)+g.dx(1,jj)*g.dxxxx(1,ii,kk,ll,mm)+...
g.dxx(1,ll,mm)*g.dxxx(1,ii,jj,kk)+g.dx(1,ll)*g.dxxxx(1,ii,jj,kk,mm)...
)*f...
+(g.dxx(1,kk,ll)*g.dxx(1,ii,jj)+...
g.dx(1,kk)*g.dxxx(1,ii,jj,ll)+...
g.dxxx(1,jj,kk,ll)*g.dx(1,ii)+...
g.dxx(1,jj,kk)*g.dxx(1,ii,ll)+...
g.dxx(1,jj,ll)*g.dxx(1,ii,kk)+...
g.dx(1,jj)*g.dxxx(1,ii,kk,ll)+...
g.dx(1,ll)*g.dxxx(1,ii,jj,kk))*obj.dx(1,mm)...
+(...
g.dxx(1,kk,mm)*g.dxx(1,ii,jj)+g.dx(1,kk)*g.dxxx(1,ii,jj,mm)+...
g.dxxx(1,jj,kk,mm)*g.dx(1,ii)+g.dxx(1,jj,kk)*g.dxx(1,ii,mm)+...
g.dxx(1,jj,mm)*g.dxx(1,ii,kk)+g.dx(1,jj)*g.dxxx(1,ii,kk,mm)...
)*obj.dx(1,ll)...
+(...
g.dx(1,kk)*g.dxx(1,ii,jj)+...
g.dxx(1,jj,kk)*g.dx(1,ii)+...
g.dx(1,jj)*g.dxx(1,ii,kk)...
)*obj.dxx(1,ll,mm)...
+(...
g.dxxx(1,jj,ll,mm)*g.dx(1,ii)+g.dxx(1,jj,ll)*g.dxx(1,ii,mm)+...
g.dxx(1,jj,mm)*g.dxx(1,ii,ll)+g.dx(1,jj)*g.dxxx(1,ii,ll,mm)...
)*obj.dx(1,kk)...
+(...
g.dxx(1,jj,ll)*g.dx(1,ii)+...
g.dx(1,jj)*g.dxx(1,ii,ll)...
)*obj.dxx(1,kk,mm)...
+g.dxx(1,jj,mm)*g.dx(1,ii)*obj.dxx(1,kk,ll)...
+g.dx(1,jj)*g.dxx(1,ii,mm)*obj.dxx(1,kk,ll)...
+g.dx(1,jj)*g.dx(1,ii)*obj.dxxx(1,kk,ll,mm)...
);
end
end

function obj=reinitialize(g)
obj=g;
obj.x=0;
oo_=g.order;
obj.dx(:)=0;
if oo_>1
obj.dxx(:)=0;
if oo_>2
obj.dxxx(:)=0;
if oo_>3
obj.dxxxx(:)=0;
if oo_>4
obj.dxxxxx(:)=0;
end
end
end
end
end

function C=autofill_symmetry(C,proto,bigset,current_level,nvars)
% autofill_symmetry --  symmetrizer for higher-order derivatives of aplanar
%
% ::
%
%
%   C=autofill_symmetry(C,proto)
%
%   C=autofill_symmetry(C,proto,bigset)
%
%   C=autofill_symmetry(C,proto,bigset,current_level)
%
%   C=autofill_symmetry(C,proto,bigset,current_level,nvars)
%
% Args:
%
%    - **C** [array]: derivatives
%
%    - **proto** [matrix]: permutations of positions of variables
%
%    - **bigset** [cell]: combination of variables
%
%    - **current_level** [integer]: level of filling of bigset
%
%    - **nvars** [integer]: number of variables
%
% Returns:
%    :
%
%    - **C** [array]: derivatives
%
% Note:
%
%    - It is assumed that the symmetric elements of the derivatives are not
%    computed. e.g. aplanar computes f(i,j,k) such that i>=j>=k. This function
%    uses that information to fill in f(i,k,j), f(j,i,k), f(j,k,i), f(k,i,j),
%    f(k,j,i), which are all equal to f(i,j,k)
%
%    - There two algorithms: a slow one, which is recursive and a fast one
%    which is non-recursive. The non-recursive algorithm is triggered when the
%    number of inputs is exactly 2
%
% Example:
%
%    See also:

siz=size(C);
order=numel(siz)-1;
if order<2
    return
end
recursive_algo=nargin>2;
if recursive_algo
    if nargin<5
        nvars=siz(2);
        if nargin<4
            current_level=1;
            if nargin<3
                proto=cell2mat(utils.gridfuncs.mypermutation(1:order));
                if nargin<2
                    bigset=cell(1,order);
                end
            end
        end
    end
    if numel(bigset)~=order
        error('number elements of bigset different from order')
    end
end

if siz(1)>1
    proto_size=siz;
    proto_size(1)=1;
    output=cell(siz(1),1);
    for irow=1:siz(1)
        if recursive_algo
            output{irow}=autofill_symmetry(reshape(C(irow,:),proto_size),...
                proto,bigset,current_level,nvars);
        else
            output{irow}=autofill_symmetry(reshape(C(irow,:),proto_size),...
                proto);
        end
    end
    C=cell2mat(output);
else
    
    if recursive_algo
        % build each set, find its variants and replace
        %----------------------------------------------
        finish=bigset{current_level};
        if isempty(finish)
            finish=nvars;
        end
        for ii=1:finish
            bigset{current_level}=ii;
            if current_level==order
                do_it();
            else
                C=autofill_symmetry(C,proto,bigset,current_level+1,nvars);
            end
        end
    else
        % find all occurrences, construct variants and replace
        %------------------------------------------------------
        [~,J,V]=find(C);
        
        [meeks{1:order+1}]=ind2sub(siz,J);
        
        nrows=numel(V);
        
        variants0=cell2mat(meeks);
        
        for irow=1:nrows
            variants=variants0(irow,2:end);
            variants=variants(proto);
            rr=size(variants,1);
            variants=mat2cell(variants,rr,ones(1,order));
            locs=unique(sub2ind(siz,ones(rr,1),variants{:}));
            C(locs)=V(irow);
        end
    end
end

    function do_it()
        bingo=C(1,bigset{:});
        if bingo
            variants0=cell2mat(bigset);
            % don't do the "diagonals" since there is only one element in
            % that position. Ideally I should remove all duplicates...
            %-------------------------------------------------------------
            if any(variants0(1)-variants0)
                variants0=variants0(proto);
                nrows=size(variants0,1);
                variants=cell(1,order+1);
                variants{1}=ones(nrows,1);
                for icol=1:order
                    variants{icol+1}=variants0(:,icol);
                end
                % do only the unique... since we do not remove the
                % duplicates above
                %----------------------------------------------------------
                locs=unique(sub2ind(size(C),variants{:}));
                C(locs)=bingo;
            end
        end
    end
end

