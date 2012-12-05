function eqtn_out=trim_symbolic_equation(eqtn_in,debug)

if nargin<2
    debug=false;
end
%%
must_change=ischar(eqtn_in);
if must_change
    eqtn_in=cellstr(eqtn_in);
end

if numel(eqtn_in)>1
    debug=false;
end
%% initialize
eqtn_out=eqtn_in;

%% white space
eqtn_out(isspace(eqtn_out))=[];

%% Vectorizable patterns and explanation
patterns={
    '(\-\-)|(\+\+)','+',                  '"--" or "++" --->"+"' 
    '(\+\-)|(\-\+)','-',                  '"+-" or "-+" ---> "-"' 
    '(?<!(\w|\)))(\+)(\w+)','$2',         '"+atom" ---> "atom" if not preceded by \w|\)' 
    '(?<!\w)(\()(\w+)(\))','$2',          '"(atom)" ---> "atom" if atom and not preceded by a function name'  
    '(\-\(\-)(\w+)(\))','$2',               '"-(-atom)" ---> "atom"'  
    '(\()(\+)(\w+)','$1$3',               '"(+atom" ---> "(atom"' 
    '(?<![\w\.])1\*','',         '"1*"--->"" if not preceeded by: \w or .'
    '[\*\^]1(?![\^/\d\.])','',   '"*1" or "^1"--->"": if not followed by \d|\^|/|\.'
    '([\d])(\.0+)(?![\d])','$1', '"digit.0"-->"digit" if not followed by \d'
    '(?<!\w)uminus\(0\)','0',    ' uminus(0) ---> 0'
    };
% (?(?=condition)(then1|then2|then3)|(else1|else2|else3))

%% now do it
if debug
    disp('original equation')
    disp(eqtn_out{1})
    old_eq=eqtn_out{1};
end
for ipat=1:size(patterns,1)
    eqtn_out=regexprep(eqtn_out,patterns{ipat,1},patterns{ipat,2});
    if debug
        disp(' ')
        disp(patterns(ipat,:))
        if isequal(old_eq,eqtn_out{1})
            disp('No effect of previous pattern')
        else
            disp(eqtn_out{1})
            old_eq=eqtn_out{1};
        end
    end
end

%% remove unncessary parentheses and zeros
for ieq=1:numel(eqtn_out)
    eqtn_out{ieq}=string_optimize.remove_unnecessary_parentheses(eqtn_out{ieq},debug);
    eqtn_out{ieq}=string_optimize.remove_zeros_and_ones(eqtn_out{ieq},debug);
end

bad=~cellfun(@isequal,eqtn_in,eqtn_out);
if any(any(bad))
    eqtn_out(bad)=trim_symbolic_equation(eqtn_out(bad),debug);
else
    if must_change
        eqtn_out=char(eqtn_out);
    end
end


end
