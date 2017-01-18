function [dec,flag]=decompose_date(x)
% annual: 1 or '1'
% bi-annual: '1990H1'
% quarterly: '1990Q1'
% monthly: '1990M1'
% weekly: '1990W1'

prototype=struct('year','','frequency','','period','','freq','');

if nargin==0
    
    dec=prototype;
    
    flag=true;
    
    return
    
end

if ischar(x) && all(all(isstrprop(x,'digit')))
    
    x=str2double(x);
    
end

if isnumeric(x)
    
    if is_serial(x)
        
        [dec]=serial2dec(x);
        
        flag=true;
        
    else
        
        [dec,flag]=decompose_yearly_date(x);
        
    end
    
elseif ischar(x)||iscellstr(x)
    
    [dec,flag]=decompose_wmqh_date(x);
    
end

end
