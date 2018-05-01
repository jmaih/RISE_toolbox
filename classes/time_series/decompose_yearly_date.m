function [dec,flag]=decompose_yearly_date(x)
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


% try turning into annual if possible and do not scream
%-------------------------------------------------------
x=char2num(x,true);

if ischar(x)||iscellstr(x)
    
    [dec,flag]=decompose_wmqh_date(x);
    
    return
    
end

if is_serial(x)
    
    [dec]=serial2dec(x,true);
    
    flag=true;
    
else
    
    dec=[];
    
    flag=false;
    
    if any(x-floor(x))
        
        return
        
        %     error('wrong annual date format')
        
    end
    
    flag=true;
    
    nx=numel(x);
    
    dec=decompose_date();
    
    dec.frequency='';
    
    dec.period=1;
    
    dec.freq=1;
    
    dec=dec(1,ones(nx,1));
    
    ii=0;
    
    while ii<nx
        
        ii=ii+1;
        
        dec(ii).year=x(ii);
        
    end
    
end

end