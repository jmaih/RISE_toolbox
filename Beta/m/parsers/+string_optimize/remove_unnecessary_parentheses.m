function eqtn=remove_unnecessary_parentheses(eqtn,debug)
if nargin<2
    debug=false;
end
% Replace "((....))" with "(....)"
% Remove unneeded parentheses
% start from the end and move back towards the beginning...
[parmatch,n]=match_parentheses(eqtn);
discard=false(1,n);
eqtn_length=length(eqtn);
for iloc=n:-1:1
    left=parmatch(1,iloc);
    right=parmatch(2,iloc);
    flagged=false;
    if (left==1 && right==eqtn_length)
        flagged=true;
    end
    if ~flagged && left-1>=1
        if any(eqtn(left-1)=='(+') && right+1<=eqtn_length
            if any(eqtn(right+1)==')+-')
                flagged=true;
            end
        end
    end
    if flagged
        discard(iloc)=true;
    end
end
discard=parmatch(:,discard);
discard=discard(:)';
if debug
    disp(upper(mfilename))
    new_eqtn=eqtn;
    new_eqtn(discard)=[];
    if isequal(eqtn,new_eqtn)
        disp('no change')
    else
        disp('final expression')
        disp(new_eqtn)
    end
end
eqtn(discard)=[];

end