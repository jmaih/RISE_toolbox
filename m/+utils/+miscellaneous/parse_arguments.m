function varargout=parse_arguments(args,varargin)
% PARSE_ARGUMENTS -- arguments parser
%
% ::
%
%
%   varargout=PARSE_ARGUMENTS(args,varargin)
%
%   struct=PARSE_ARGUMENTS(args,struct)
%
%   struct=PARSE_ARGUMENTS(args,struct,'-hopover')
%
% Args:
%
%    - **args** [n x 3 or n x 4 cell]:
%      - column 1: arguments' names
%      - column 2: defaults values for the arguments
%      - column 3: function handle checking the value of the arguments
%      - column 4: message to be issued in case of error
%
%    - **varargin** []: pairwise arguments. Each first argument is a char and
%    the second is a value to be parsed for correctness
%
%    - **struct** [struct]: structure with the options to parse
%
%    - **'-hopover'** []: if present, the options that are not found in args
%    are skipped. Else an error is issued
%
% Returns:
%    :
%
%    - **varargout** [list|struct]: if the second input argument is a struct,
%    a struct is return in varargout. Else a list is returned. In both cases,
%    only the elements relating to the function of interest (i.e. the elements
%    in **args**) are returned.
%
% Note:
%
% Example:
%
%    See also:

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

hopover=~isempty(varargin) &&...
    ischar(varargin{end})&& ...
    strcmp(deblank(varargin{end}),'-hopover');

if hopover
    
    varargin(end)=[];
    
end

struct_style=length(varargin)==1 && isstruct(varargin{1});

if hopover && ~struct_style
    
    error('hopover can only be used with structures')
    
end

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
    
    if ~any(prop_loc)
        
        if hopover
            
            continue
            
        else
            
            error(sprintf('%s is not a valid property',prop_names{iprop})) %#ok<SPERR>
            
        end
        
    end
    
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