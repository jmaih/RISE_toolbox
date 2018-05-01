function Value= get(obj,varargin)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:


%  get  Access/Query time series property values.
%
%     VALUE = get(TS,'PropertyName') returns the value of the
%     specified property of the time series object.  An equivalent
%     syntax is
%
%         VALUE = TS.PropertyName
%
%     get(TS) displays all properties of TS and their values.

if isempty(obj)
    
    Value = [];
    
    return
    
end

if numel(obj)==1 % get on a scalar timeseries
    
    Value = uttsget(obj,varargin{:});
    
    return
    
end

% Process array values
if nargin>=2 % get on a timeseries array with specified properties
    
    if ischar(varargin{1})
        
        Value = cell(size(obj));
        
        for k=1:numel(obj)
            
            Value{k} = uttsget(obj(k),varargin{:});
        
        end
        
    elseif iscell(varargin{1})
        
        props = varargin{1};
        
        Value = cell(numel(obj),length(props));
        
        for k=1:numel(obj)
            
            for j=1:length(props)
                
                Value{k,j} = uttsget(obj(k),props{j});
            
            end
            
        end
        
    end
    
else % Return a stuct array for a timeseries array with no props
    
    for k=numel(obj):-1:1
        
        Value(k) = uttsget(obj(k),varargin{:});
    
    end
    
end

end

function ValueOut = uttsget(h,varargin)
%
% ts utility function

ni = nargin;

narginchk(1,2);

if ni==2,
    % GET(H,'Property') or GET(H,{'Prop1','Prop2',...})
    Property = varargin{1};
    
    CharProp = ischar(Property);
    
    if CharProp,
        
        Property = {Property};
    
    elseif ~iscellstr(Property)
        
        error('invalid property')
    
    end
    
    % Loop over each queried property
    Nq = numel(Property);
    
    Value = cell(length(h),Nq);
    
    for i=1:Nq
        % Find match for k-th property name and get corresponding value
        % RE: a) Must include all properties to detect multiple hits
        %     b) Limit comparison to first 7 chars (because of iodelaymatrix)
        try
            
            if numel(h)==1 % Do not index into timeseries - they are not really arrays
                
                Value{1,i} = h.(Property{i});
            
            else
                
                for k=1:numel(h)
                    
                    Value{k,i} = h(k).(Property{i});
                
                end
                
            end
            
        catch me
            
            rethrow(me)
        
        end
        
    end
    
    % Strip cell header if PROPERTY was a string.
    if CharProp,
        
        ValueOut = Value{1};
    
    else
        
        ValueOut = Value;
    
    end
    
else
    
    if all(~ishandle(h)) % This is an MCOS object/array
        
        classH = metaclass(h);
        
        Public=strcmp({classH.PropertyList.GetAccess},'public');
        
        Hidden=[classH.PropertyList.Hidden];
        
        good=Public & ~Hidden;
        
        PropNames = {classH.PropertyList.Name};
        
        PropNames=PropNames(good);
        
    else
        
        classH = classhandle(h(1));
        
        propH  = classH.Properties(:);
        
        PropNames = {};
        
        for k=1:length(propH)
            
            if strcmp(propH(k).AccessFlags.PublicGet,'on') && strcmp(propH(k).Visible,'on')
                
                PropNames = [PropNames; {propH(k).Name}]; %#ok<AGROW>
                
            end
            
        end
        
    end
    
    if numel(h)>1
        
        for j=numel(h):-1:1
            
            for k=length(PropNames):-1:1
                
                PropValues{k} = h(j).(PropNames{k});
                
            end
            
            Value(j) = cell2struct(PropValues,PropNames,2);
            
        end
        
        ValueOut = Value;
        
    else
        
        for k=length(PropNames):-1:1
            
            PropValues{k} = h.(PropNames{k});
            
        end
        
        ValueOut = cell2struct(PropValues,PropNames,2);
        
    end
    
end

end