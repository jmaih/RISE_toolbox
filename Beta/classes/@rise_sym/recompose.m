function c=recompose(func,varargin)
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