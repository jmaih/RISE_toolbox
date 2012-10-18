function a=tryevaluate(a)
% checks whether a string can be evaluated
allowed_set='0123456789.+-*()/^';

method=3;
switch method
    case 1
        if all(ismember(a,allowed_set))%<--- ~any(isletter(a))
            a=num2str(eval(a),10);
        end
    case 2
        flag=true;
        na=length(a);
        for ii=1:length(a)
            flag=flag && any(a(ii)==allowed_set);
            if ~flag
                break
            end
        end
        if flag && na>3 % if it is already a 
            a=num2str(eval(a),10);
        end
    case 3
        flag=~any(isstrprop(a,'alpha'));
        if flag
            cntrl=a;
            cntrl(isstrprop(cntrl,'digit'))=[];
            flag=~isempty(cntrl) && ~isequal(cntrl,'.');
            if flag
                flag=false;
                for ii=1:length(cntrl)
                    if any(cntrl(ii)=='+-*^/')
                        flag=true;
                        break
                    end
                end
                if flag
                    a=num2str(eval(a),10);
                end
            end
        end
    otherwise
        error([mfilename,':: case not implemented'])
end
end

