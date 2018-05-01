function d=splice(d,varargin)

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

if isa(d,'ts')
    
    d=pages2struct(d);
    
end

check_struct(d)

if length(varargin) > 1
    
    for i = 1 : length(varargin)
        
        d = splice(d,varargin{i});
        
    end
    
else
    
    do_it()
    
end

%--------------------------------------------------------------------------

    function do_it()
        
        s = varargin{1};
        
        if isa(s,'ts')
            
            s=pages2struct(s);
            
        end
        
        check_struct(s)
        
        dnames = fieldnames(d);
        
        snames = fieldnames(s);
        
        all_names = union(dnames,snames);
        
        for j = 1 : numel(all_names)
            
            if ~isfield(s,all_names{j})
                
                continue
                
            end
            
            if ~isfield(d,all_names{j})
                
                d.(all_names{j}) = s.(all_names{j});
                
                continue
                
            end
            
            x = d.(all_names{j});
            
            y = s.(all_names{j});
            
            if isa(x,'ts') && isa(y,'ts')
                
                if isequal(x.frequency,y.frequency)
                    
                    % Two non-empty ts with the same frequency.
                    d.(all_names{j}) = tsvercat(x,y);
                    
                elseif isempty(x.data)
                    
                    % Two empty ts or the first non-empty and the
                    % second empty; use the first input anyway.
                    d.(all_names{j}) = y;
                    
                elseif isempty(y.data)
                    
                    % Only the second ts is non-empty.
                    d.(all_names{j}) = x;
                    
                else
                    
                    % Two non-empty ts with different frequencies.
                    d.(all_names{j}) = x;
                    
                end
                
            else
                
                % At least one non-ts input, use the second input.
                d.(all_names{j}) = y;
                
            end
            
        end
        
        tempList = fieldnames(s);
        
        tempList = setdiff(tempList,all_names);
        
        for j = 1 : length(tempList)
            
            d.(tempList{j}) = s.(tempList{j});
            
        end
        
    end

end

function check_struct(x)

if ~isstruct(x) % || any(~cellfun(@isstruct,varargin))
    
    error('All input arguments must be structure of time series (ts) or time series');
    
end

end

function v=tsvercat(varargin)

nargs=length(varargin);

tmp=struct();

for iarg=1:nargs
    
    x=varargin{iarg};
    
    dx=x.date_numbers;
    
    tmp(iarg).start=dx(1);
    
    tmp(iarg).end=dx(end);
    
    tmp(iarg).data=x.data;
    
    tmp(iarg).r=size(x.data,1);
    
    tmp(iarg).c=size(x.data,2);
    
    tmp(iarg).p=size(x.data,3);
    
end

s=min([tmp.start]);

b=max([tmp.end]);

c=max([tmp.c]);

p=max([tmp.p]);

d=s:b;

v=nan(numel(d),c,p);

for iarg=1:nargs
    
    fill_it()
    
end

v=ts(d(1),v);

    function fill_it()
        
        startx=find(tmp(iarg).start==d);
        
        if tmp(iarg).c<c && tmp(iarg).c==1
            
            tmp(iarg).data=tmp(iarg).data(:,ones(1,c));
            
            tmp(iarg).c=c;
            
        end
        
        v(startx:startx+tmp(iarg).r-1,1:tmp(iarg).c,1:tmp(iarg).p)=...
            tmp(iarg).data;
        
    end

end