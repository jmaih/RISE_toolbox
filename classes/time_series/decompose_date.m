function [dec,flag]=decompose_date(x)
% annual: 1 or '1'
% bi-annual: '1990H1'
% quarterly: '1990Q1'
% monthly: '1990M1'
% weekly: '1990W1'

prototype=struct('year','','WMQH','','period','','frequency','','freq','');

if nargin==0
    
    dec=prototype;
    
    flag=true;
    
    return
    
end


if ischar(x)
    
    x=cellstr(x);
    
end

if isnumeric(x)
    
    % serial numbers are treated in a different way
    %----------------------------------------------
    
    ss=is_serial(x);
    
    if all(ss)
        
        x=serial2date(x);
        
        [dec,flag]=decompose_date(x);
        
        return
        
    elseif any(ss)
        
        error('mixing serial numbers with dates not allowed')
        
    end
    
    x0=x;
    
    x=cell(size(x0));
    
    for iii=1:numel(x0)
                
        x{iii}=num2str(x0(iii));
        
    end
    
    
end

if ~iscellstr(x)
    
    error('inputs must be numeric, char or cellstr')
    
end

dec=[];

% row vector
old_size=size(x);

x=upper(x(:).');

nx=numel(x);

tmp=prototype(1,ones(nx,1));

% separate annual dates from the rest
%------------------------------------
is_annual=cellfun(@(x)all(isstrprop(x,'digit')),x,'uniformOutput',true);

success=true;

if any(is_annual)
    
    [do_this,flag]=decompose_yearly_date(x(is_annual));
    
    if flag
        
        tmp(is_annual)=do_this;
        
    else
        
        success=false;
        
    end
    
end
    
if success && any(~is_annual)
    
    [do_this,flag]=decompose_wmqh_date(x(~is_annual));
    
    if flag
        
        tmp(~is_annual)=do_this;
        
    else
        
        success=false;
        
    end
    
end

if success
    
    dec=reshape(populate_extra_fields(tmp),old_size);
    
end

end

function tmp=populate_extra_fields(tmp)

fmap=frequency_map();

quick_struct=struct();

for istr=1:numel(fmap.strings)-1 % don't take the year
    
    quick_struct.(fmap.strings{istr})=fmap.code(istr);
    
end

for icase=1:numel(tmp)
    
    frequency=tmp(icase).WMQH;
    
    tmp(icase).frequency=frequency;
    
    if isempty(frequency)
        
        tmp(icase).freq=1;
        
    else
        
        tmp(icase).freq=quick_struct.(frequency);
        
    end
    
end

end
