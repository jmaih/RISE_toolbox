function flag=is_date(x)
% annual: 1 or '1'
% bi-annual: '1990H1'
% quarterly: '1990Q1'
% monthly: '1990M1'
% weekly: '1990W1'

if ischar(x)
    
    x=cellstr(x);
    
end

if isnumeric(x)
    
    x0=x;
    
    x=cell(size(x0));
    
    for iii=1:numel(x0)
        
        x{iii}=num2str(x0(iii));
        
    end
    
end

if numel(x)>1
    
    flag=false(size(x));
    
    for iii=1:numel(x)
        
        flag(iii)=is_date(x(iii));
        
    end
    
    return
    
end

flag=false;

if ~iscellstr(x)
    
    error('inputs must be numeric, char or cellstr')
    
end

x=upper(x);

% separate annual dates from the rest
%------------------------------------
is_annual=cellfun(@(x)all(isstrprop(x,'digit')),x,'uniformOutput',true);

n_ann=sum(is_annual);

if n_ann
    
    % this should be automatically true
    %----------------------------------
    
    [~,flag]=decompose_yearly_date(x(is_annual));
    
    if ~flag||numel(is_annual)==n_ann
        
        return
        
    end
    
end

if any(~is_annual)
    
    [~,flag]=decompose_wmqh_date(x(~is_annual));
    
end

end
