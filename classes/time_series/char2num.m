function [dat,is_colon]=char2num(dat,silent)
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


is_colon=false;

if isnumeric(dat)
    
    return
    
end

if ischar(dat)
    
    colon_pos=find(dat==':');
    
    if numel(colon_pos)>1
        
        error('Only a max of one ":" is allowed in dates')
        
    end
    
    is_colon=~isempty(colon_pos);
    
end

trysplit=dat;

if ~iscellstr(trysplit)
    
    if is_colon
        
        trysplit=regexp(trysplit,':','split');
        
    else
        
        trysplit={trysplit};
        
    end
    
end

is_annual=cellfun(@(x)all(isstrprop(x,'digit')),trysplit,...
    'uniformOutput',true);

if is_annual
    
    dat=cellfun(@(x)str2double(x),trysplit,'uniformOutput',true);
    
else
    
    if ~silent
        
        error('I could not turn char to double')
        
    end
    
end
end