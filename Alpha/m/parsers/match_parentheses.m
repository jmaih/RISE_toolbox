function [parmatch,n]=match_parentheses(eqtn,type)
if nargin<2
    type='parentheses';
end
switch lower(type)
    case {'parent','parenthesis','parentheses'}
        left_par='(';
        right_par=')';
    case {'bracket','brackets'}
        left_par='[';
        right_par=']';
    case {'braces','brace','curly_braces','curly_brackets','curlybraces','curlybrackets'}
        left_par='{';
        right_par='}';
    otherwise
        error(['unknown type ',type,'. second argument expected to be ''parentheses'', ''brackets'' or ''braces'''])
end

test=true;
if test
    str_length=length(eqtn);
    iter=0;
    locs=false(1,str_length);
    while iter<str_length
        iter=iter+1;
        if eqtn(iter)==left_par||eqtn(iter)==right_par
            locs(iter)=true;
        end
    end
    locs=find(locs);
else
    locs=regexp(eqtn,['\',left_par,'|\',right_par],'start');
end

left_right=eqtn(locs);
n=0.5*numel(left_right);
if n~=floor(n)
    error('unbalanced number of braces, brackets or parentheses')
end
parmatch=nan(2,n);
opened_nbr=0;
closed_nbr=0;
opened_list=[];
for ii=1:numel(locs)
    if left_right(ii)==left_par %<--- strcmp(left_right(ii),left_par)
        opened_nbr=opened_nbr+1;
        parmatch(1,opened_nbr)=locs(ii);
        opened_list=[opened_list,opened_nbr]; %#ok<AGROW>
    else
        parmatch(2,opened_list(end))=locs(ii);
        opened_list(end)=[];
        closed_nbr=closed_nbr+1;
    end
    if closed_nbr>opened_nbr
        error(['closed ',left_right(ii),' at position ',int2str(locs(ii)),' was not opened'])
    end
end
if ~isempty(opened_list)
    disp(locs(opened_list))
    error('open parenthesis, bracket or brace in positions listed were not closed')
end

% OLD and slow algorithm
% left_pars=find(eqtn=='(');
% n=numel(left_pars);
% parmatch=nan(2,n);
% if n
%     right_pars=find(eqtn==')');
%     if numel(right_pars)~=n
%         error('number of left parentheses does not match')
%     end
%     for ip=n:-1:1
%         loc=find(right_pars>left_pars(end),1,'first');
%         if isempty(loc)
%             error('parentheses mismatch')
%         end
%         parmatch(1,ip)=left_pars(end);
%         parmatch(2,ip)=right_pars(loc);
%         left_pars(end)=[];
%         right_pars(loc)=[];
%     end
% end
end

