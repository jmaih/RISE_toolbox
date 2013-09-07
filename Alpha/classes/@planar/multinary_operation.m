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
    if isa(varargin{ivar},'planar')
        numargs{ivar}=varargin{ivar}.func;
    else
        numargs{ivar}=varargin{ivar};
    end
end
if numflag
    obj=feval(fcn,numargs{:});
else
    obj=planar(fcn,varargin);
    obj=set(obj,'incidence',get_incidence(varargin{:}));
%     obj=commit(obj);
end
end

function occur=get_incidence(varargin)
occur=[];
for iv=1:length(varargin)
    if isa(varargin{iv},'planar') && ~isempty(varargin{iv}.incidence)
        if isempty(occur)
            occur=varargin{iv}.incidence;
        else
            occur=occur|varargin{iv}.incidence;
        end
    end
end
end