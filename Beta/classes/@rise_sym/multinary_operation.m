function obj=multinary_operation(fcn,varargin)
numflag=true;
nargs=length(varargin);
for ivar=1:nargs
    numflag=numflag && isnumeric(varargin{ivar});
    if ~numflag
        break
    end
    if ivar==1
        numargs=cell(1,nargs);
    end
    if isa(varargin{ivar},'rise_sym')
        numargs{ivar}=varargin{ivar}.func;
    else
        numargs{ivar}=varargin{ivar};
    end
end
if numflag
    obj=feval(fcn,numargs{:});
else
    obj=rise_sym(fcn,varargin,rise_sym.get_incidence(varargin{:}));
    obj=commit(obj);
end
end
%{
function obj=unary_operation(fcn,a)
if isnumeric(a.func)
    obj=feval(fcn,a);
else
    obj=rise_sym(fcn,{a},a.incidence);
    obj=commit(obj);
end
end
function obj=binary_operation(fcn,a,b)
if isnumeric(a.func) && isnumeric(b.func)
    obj=feval(fcn,a,b);
else
    obj=rise_sym(fcn,{a,b},rise_sym.get_incidence(a,b));
    obj=commit(obj);
end
end
function obj=trinary_operation(fcn,a,b,c)
if isnumeric(a.func) && isnumeric(b.func) && isnumeric(c.func)
    obj=feval(fcn,a,b,c);
else
    obj=rise_sym(fcn,{a,b,c},rise_sym.get_incidence(a,b,c));
    obj=commit(obj);
end
end
%}
