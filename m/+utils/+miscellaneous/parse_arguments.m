function varargout=parse_arguments(args,varargin)
% parse_arguments -- arguments parser
%
% Syntax
% -------
% ::
%
%   varargout=parse_arguments(args,varargin)
%
%   struct=parse_arguments(args,struct)
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

% {names,defaults,check_funcs,error_msg}
assert(iscell(args) && any(size(args,2)==[3,4]),...
    'the first input argument must be a cell array with 3 or 4 columns')
arg_names=args(:,1).';
defaults=args(:,2).';
checks=args(:,3).';
error_msg=args(:,1:0).';
if size(args,2)==4
    error_msg=args(:,end).';
    assert(iscellstr(error_msg),'fourth column of first input must be a cellstr')
end
struct_style=length(varargin)==1 && isstruct(varargin{1});
nout=nargout;
nargs=numel(arg_names);
if struct_style
    if nout~=1
        error('when the second input is a struct, the output must be a struct')
    end
    prop_names=fieldnames(varargin{1});
    prop_vals=struct2cell(varargin{1});
    varargout={cell2struct(defaults,arg_names,2)};
else
    if nout>nargs
        error('number of output arguments cannot exceed the number of rows of the first input')
    end
    varargout=defaults;
    assert(rem(length(varargin),2)==0,'arguments must come in pairs')
    prop_names=varargin(1:2:end);
    prop_vals=varargin(2:2:end);
end
assert(iscellstr(arg_names),'first column of first input must be a cellstr')
for iprop=1:numel(prop_names)
    prop_loc=strcmp(prop_names{iprop},arg_names);
    assert(any(prop_loc),...
        sprintf('%s is not a valid property',prop_names{iprop}))
    value=prop_vals{iprop};
    if isempty(value) && ~struct_style
        continue
    end
    if isempty(error_msg)
        assert_args={checks{prop_loc}(value)};
    else
        assert_args={checks{prop_loc}(value),error_msg{prop_loc}};
    end
    assert(assert_args{:})
    if struct_style
        varargout{1}.(prop_names{iprop})=value;
    else
        varargout{prop_loc}=value;
    end
end

end