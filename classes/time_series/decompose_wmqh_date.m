function [dec,flag]=decompose_wmqh_date(x)
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

if ~(ischar(x)||iscellstr(x))
    
    error('input must be a char or a cellstr')
    
end
% try turning into annual if possible and do not scream
%-------------------------------------------------------
x=char2num(x,true);

if isnumeric(x)
    
    [dec]=serial2dec(x);
    
    return
    
end

flag=false;

dec=[];

fmap=frequency_map();

WMQH=strrep(cell2mat(strcat(fmap.strings','|')),'||','');

test=regexpi(x,['\d+?(',WMQH,')\d+|\d+'],'tokens');

test=[test{:}];

if iscellstr(x) % iscell(test{1})
    
    test=[test{:}];
    
end

% first test
%-----------
nonsense=cellfun(@isempty,test,'uniformOutput',true);

if any(nonsense)||isempty(nonsense)
    
    return
    
end

express='(?<year>\d+)(?<frequency>(W|M|Q|H))(?<period>\d+)';

tmp=regexpi(x,express,'names');

if iscellstr(x) % iscell(test{1})
    
    tmp=[tmp{:}];
    
end

for id=1:numel(tmp)
    
    tmp(id).year=str2double(tmp(id).year);
    
    tmp(id).period=str2double(tmp(id).period);
    
    i_check=is_valid_wmqh_periodicity();
    
    if ~i_check
        
        break
        
    end
    
end

if i_check
    
    freq=[tmp.freq];
    
    if all(freq==freq(1))
        
        flag=true;
        
        dec=tmp;
        
    end
    
end

    function i_check=is_valid_wmqh_periodicity()
        
        per=tmp(id).period;
        
        i_check=per>0 && per==floor(per);
        
        if ~i_check
            
            return
            
        end
        
        switch upper(tmp(id).frequency)
            
            case 'W'
                
                i_check=per<=52;
                
                tmp(id).freq=52;
                
            case 'M'
                
                i_check=per<=12;
                
                tmp(id).freq=12;
                
            case 'Q'
                
                i_check=per<=4;
                
                tmp(id).freq=4;
                
            case 'H'
                
                i_check=per<=2;
                
                
                tmp(id).freq=2;
        end
        
        
    end


end
