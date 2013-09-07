function [separators,len,nsep,delimiters]=find_separators(eqtn,delimiters)
if nargin<2
    delimiters='';
end
if isempty(delimiters)
    delimiters='/*-+^(){}[]';
end
if ~ischar(eqtn)
    error('eqtn must be a char object')
end
len=length(eqtn);
separators= cell2mat(strcat('\',num2cell(delimiters)));
% location of delimiters
separators= regexp(eqtn,['[',separators,']'],'start');
seps=eqtn(separators);
nsep=numel(seps);
separators=[
    nan(1,nsep)
    separators
    nan(1,nsep)
    false(1,nsep)
    ];
for idel=1:length(delimiters)
    % type of delimiter
    separators(1,seps==delimiters(idel))=idel;
end
open_pars=[];
open_braces=[];
open_brackets=[];
for jsep=1:length(seps)
    if separators(1,jsep)==6
        open_pars=[open_pars,separators(2,jsep)]; %#ok<*AGROW>
        if separators(2,jsep)>1 && any(eqtn(separators(2,jsep)-1)==['a':'z','_','A':'Z','0':'9'])
            % parenthesis of function
            separators(4,jsep)=true;
        end
    elseif separators(1,jsep)==7
        if isempty(open_pars)
            error(['closing parenthesis not opened at position ',int2str(iter)])
        end
        last_sep=open_pars(end);
        % closing of parenthesis
        separators(3,separators(2,:)==last_sep)=separators(2,jsep);
        open_pars(end)=[];
    elseif separators(1,jsep)==8
        open_braces=[open_braces,separators(2,jsep)]; %#ok<*AGROW>
        if separators(2,jsep)>1 && any(eqtn(separators(2,jsep)-1)==['a':'z','_','A':'Z','0':'9'])
            % parenthesis of function
            separators(4,jsep)=true;
        end
    elseif separators(1,jsep)==9
        if isempty(open_braces)
            error(['closing brace not opened at position ',int2str(iter)])
        end
        last_sep=open_braces(end);
        % closing of parenthesis
        separators(3,separators(2,:)==last_sep)=separators(2,jsep);
        open_braces(end)=[];
    elseif separators(1,jsep)==10
        if separators(2,jsep)>1 && any(eqtn(separators(2,jsep)-1)==['a':'z','_','A':'Z','0':'9'])
            error(['brackets not allowed as functions. check position ',int2str(iter)])
        end
        open_brackets=[open_brackets,separators(2,jsep)]; %#ok<*AGROW>
    elseif separators(1,jsep)==11
        if isempty(open_brackets)
            error(['closing brace not opened at position ',int2str(iter)])
        end
        last_sep=open_brackets(end);
        % closing of parenthesis
        separators(3,separators(2,:)==last_sep)=separators(2,jsep);
        open_brackets(end)=[];
    end
end
openings=separators(1,:)==6|separators(1,:)==8|separators(1,:)==10;
sepOpen=separators(:,openings);
not_closed=isnan(sepOpen(3,:));
if any(not_closed)
    disp(sepOpen(2,not_closed))
    error('Parentheses, braces or brackets not closed')
end
end
