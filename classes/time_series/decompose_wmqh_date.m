function [dec,flag]=decompose_wmqh_date(x)

if ~iscellstr(x)
    
    error('input must be a cellstr')
    
end

dec=[];

fmap=frequency_map();

WMQH=strrep(cell2mat(strcat(fmap.strings','|')),'||','');

test=regexp(x,['\d+?(',WMQH,')\d+|\d+'],'tokens');

test=[test{:}];

% first test
%-----------
nonsense=cellfun(@isempty,test,'uniformOutput',true);

if any(nonsense)
    
    flag=false;
    
else
    
    test=[test{:}];
    
    express=['(?<year>\d+)(?<WMQH>',...
        parser.cell2matize(unique(test)),...
        ')(?<period>\d+)'];
    
    tmp=regexp(x,express,'names');
    
    tmp=[tmp{:}];
    
    for id=1:numel(tmp)
        
        i_check=is_valid_wmqh_periodicity(tmp(id));
        
        if ~i_check
            
            flag=false;
            
            break
            
        end
        
    end
    
    if i_check
        
        flag=true;
        
        dec=reshape(add_remaining_fields(tmp),size(x));
        
    end
    
end

end

function dec=add_remaining_fields(dec)

tmp_=decompose_date();

tmp_=tmp_(ones(size(dec)));

subset=fieldnames(dec);

for iset=1:numel(subset)
    
    name=subset{iset};
    
    for icase=1:numel(dec)
        
        tmp_(icase).(name)=dec(icase).(name);
        
    end
    
end

dec=tmp_;

end

function i_check=is_valid_wmqh_periodicity(tmp)

per=str2double(tmp.period);

i_check=per>0;

if ~i_check
    
    return
    
end

switch tmp.WMQH
    case 'W'
        i_check=per<=52;
    case 'M'
        i_check=per<=12;
    case 'Q'
        i_check=per<=4;
    case 'H'
        i_check=per<=2;
end

end
