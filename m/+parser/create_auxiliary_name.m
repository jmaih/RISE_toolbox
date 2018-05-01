function new_name=create_auxiliary_name(name,lead_or_lag,add_prefix)
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

if nargin<3
    
    add_prefix=false;
    
end

if ~(ischar(name)||iscellstr(name))
    
    error('first argument must be a char or a cellstr')

end

if abs(lead_or_lag)>0
    
    if sign(lead_or_lag)==1
        
        type='LEAD';
    
    elseif sign(lead_or_lag)==-1
        
        type='LAG';
    
    end
    
    new_name=sprintf('%s_%0.0f_%s',type,abs(lead_or_lag),name);

else
    
    if add_prefix
        
        new_name=sprintf('AUX_%s',name);
        
    else
        
        new_name=name;
        
    end
    
end

end