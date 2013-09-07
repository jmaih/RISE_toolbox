function [eqtn,obj]=recompose(obj)
Args=obj.args;
for iarg=1:numel(Args)
    if isa(obj.args{iarg},'rise_sym')
        obj.args{iarg}.ncalls=obj.args{iarg}.ncalls+1;
        % the effective number of calls might be lower if the argument does
        % not survive. But this is some kind of upper bound
        % handle will push this information automatically.
    else
        obj.args{iarg}=rise_sym(obj.args{iarg});
    end
    Args{iarg}=get_reference();
end
eqtn=compact_form(obj.func,Args{:});

    function [c]=get_reference()% ,obj
        if ~isempty(obj.args{iarg}.ref)
            c=obj.args{iarg}.ref;
        elseif isnumeric(obj.args{iarg}.func)
            c=sprintf('%0.16g',obj.args{iarg}.func);
        elseif isempty(obj.args{iarg}.args)
            c=obj.args{iarg}.func;
        else
            disp([mfilename,':: this case should never happen, please contact junior.maih@gmail.com'])
            obj.args{iarg}=commit(obj.args{iarg});
            c=obj.args{iarg}.ref;
        end
        
    end
end

function c=compact_form(func,varargin)
switch func
    % unary functions
    case 'uminus'
        c=['-',varargin{1}];
    case {'sqrt','abs','log','log10','exp','cos','sin','tan','acos',...
            'asin','atan','cosh','sinh','tanh','acosh','asinh','atanh',...
            'sign','erf'}
        c=[func,'(',varargin{1},')'];
        % binary functions
    case 'plus'
        c=[varargin{1},'+',varargin{2}];
    case 'minus'
        c=[varargin{1},'-',varargin{2}];
    case 'mtimes'
        c=[varargin{1},'*',varargin{2}];
    case 'mrdivide'
        c=[varargin{1},'/',varargin{2}];
    case 'mpower'
        c=[varargin{1},'^',varargin{2}];
    case {'max','min'}
        c=[func,'(',varargin{1},',',varargin{2},')'];
    case 'gt'
        c=[varargin{1},'>',varargin{2}];
    case 'ge'
        c=[varargin{1},'>=',varargin{2}];
    case 'eq'
        c=[varargin{1},'==',varargin{2}];
    case 'lt'
        c=[varargin{1},'<',varargin{2}];
    case 'le'
        c=[varargin{1},'<=',varargin{2}];
    case 'ne'
        c=[varargin{1},'~=',varargin{2}];
    case 'and'
        c=[varargin{1},'&',varargin{2}];
    case 'or'
        c=[varargin{1},'|',varargin{2}];
        % trinary functions
    case {'normalpdf','normalcdf'}
        c=[func,'(',varargin{1},',',varargin{2},',',varargin{3},')'];
    case 'if_then_else'
        if strcmp(varargin{3},'0')
            c=[varargin{1},'*',varargin{2}];
        else
            c=[varargin{1},'*',varargin{2},'+(1-',varargin{1},')*',varargin{3}];
        end
    case 'if_elseif'
        tmp=cell2mat(strcat(varargin(:)',','));
        tmp(isspace(tmp))=[];
        c=['if_elseif(',tmp(1:end-1),')'];
    otherwise
        error(['unrecognized function ',func])
end
end