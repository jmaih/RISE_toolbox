function Out = set(this,varargin)
%SET  Set properties of time series object.
%
%   SET(THIS,'PropertyName',VALUE) sets the property 'PropertyName'
%   of the time series THIS to the value VALUE.  An equivalent syntax
%   is
%       THIS.PropertyName = VALUE
%
%   SET(THIS,'Property1',Value1,'Property2',Value2,...) sets multiple
%   time series property values with a single statement.
%
%   SET(THIS,'Property') displays values for the specified property in THIS.
%
%   SET(THIS) displays all properties of THIS and their values.
%
%   See also TS\GET.

% N.B: Add a no_check property to skip the checking !!!!!!
%---------------------------------------------------------

ni = nargin;

no = nargout;

if builtin('isempty',this)
    
    this = ts;
    
end

c = metaclass(this);

AllProps = {c.PropertyList.Name};%properties(this);

bad=[c.PropertyList.Dependent]|[c.PropertyList.Constant]|[c.PropertyList.Hidden];

% Get public properties and their assignable values
AllProps=AllProps(~bad);

if ni<=2
    
    PropValues = cell(length(AllProps),1);
    
    for k=1:length(AllProps)
        
        PropValues{k} = this.(AllProps{k});
        
    end
    
end

% Handle read-only cases
if ni==1
    % SET(THIS) or S = SET(THIS)
    if numel(this)~=1
        
        error('when the number of input is 1, only one object should be set');
        
    end
    
    if no
        
        Out = cell2struct(PropValues,AllProps,1);
        
    else
        
        disp(cell2struct(PropValues,AllProps,1))
        
    end
    
elseif ni==2
    % SET(THIS,'Property') or STR = SET(THIS,'Property')
    % Return admissible property value(s)
    if numel(this)~=1
        
        error('Only one object may be passed to set for information about possible property values.');
        
    end
    
    imatch = strcmpi(varargin{1},AllProps);
    
    if no
        
        Out = PropValues{imatch};
        
    else
        
        disp(PropValues{imatch})
        
    end
    
else
    % SET(THIS,'Prop1',Value1, ...)
    if rem(ni-1,2)~=0
        
        error('Property/value pairs must come in even number.')
        
    end
    
    for i=1:2:ni-1
        
        propName = tspnmatch(varargin{i},AllProps);
        
        % Validate that base @ts property values are not structs
        if isstruct(varargin{i+1})
            
            for k=1:length(c.PropertyList)
                
                if strcmpi(c.PropertyList(k).Name,propName) && ...
                        isequal(c.PropertyList(k).DefiningClass,?ts)
                    
                    error('Time series properties cannot be assigned to structures')
                    
                end
                
            end
            
        end
        
        for k=1:numel(this)
            
            this(k).(char(propName)) = varargin{i+1};
            
        end
        
    end
    
    % Assign this in caller's workspace
    tsname = inputname(1);
    
    this=check_consistency(this);
    
    if no
        
        Out = this;
        
    elseif ~isempty(tsname)
        
        assignin('caller',tsname,this)
        
    else
        
        warning(['Cannot assign properties in-place when a timeseries is ',...
            'defined using an expression or a sub-referenced array.']);
        
    end
    
end

end


function [Property,imatch] = tspnmatch(Name,PropList)

if (~ischar(Name) || size(Name,1)>1) && ~(isstring(Name) && isscalar(Name))
    
    error('Wrong specification of property name')
    
end

% Find all matches
imatch = find(strcmpi(Name,PropList));

% Get matching property name
switch length(imatch)
    
    case 1
        % Single hit
        Property = PropList{imatch};
        
    case 0
        % No hit
        error(['property with name "', Name,'" not found or not settable'])
        
    otherwise
        
        imatch = imatch(1);
        
end

end


function this=check_consistency(this)

siz=ts.check_size(this.data);

nvars=siz(2);

% annual dates
%-------------
if isnumeric(this.start)
    
    if is_serial(this.start)
        
        this.start=serial2date(this.start);
        
    elseif ceil(this.start)==floor(this.start) % annual dates
        
        this.start=int2str(this.start);
        
    end
    
end

this.varnames=check_conformity(this.varnames,'varnames');

this.description=check_conformity(this.description,'description');

    function a=check_conformity(a,type_)
        
        na=numel(a);
        
        if all(cellfun(@isempty,a,'uniformOutput',true))
            
            if na~=nvars
                
                a=repmat({''},1,nvars);
                
            end
            
        elseif na~=nvars
            
            error(['mismatch between # columns and #',type_])
            
        end
        
    end

end