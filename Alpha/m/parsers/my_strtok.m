function [token,remainder,start,finish]=my_strtok(eqtn,delimiters)
start=[];
finish=[];
ii=1;
len = length(eqtn);
token = ''; remainder = '';
goback=false;
while any(eqtn(ii) == delimiters)
    ii = ii + 1;
    if ii > len
        goback=true;
        break
    end
end
if ~goback
    start = ii;
    while ~any(eqtn(ii) == delimiters)
        ii = ii + 1;
        if ii > len
            break
        end
    end
    finish = ii - 1;
    token = eqtn(start:finish);
    remainder = eqtn(finish + 1:len);
end
%     delims=eqtn(1:start-1)
end


