function ss=char2serial(x)
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

if iscellstr(x)
    if numel(x)>1
        error('I cannot process more than one element at a time')
    end
    x=char(x);
end
if size(x,1)>1
    error('I cannot process more than one element at a time')
end
x(isspace(x))=[];
ss=[];
while ~isempty(x)
    [tok,rest]=strtok(x,':,');
    if ~isempty(tok)
        stok=date2serial(tok);
        if strcmp(x(1),':')
            if isempty(ss)
                error('wrong date specification')
            end
            freq_old=serial2frequency(ss(end));
            freq_new=serial2frequency(stok);
            if freq_old~=freq_new
                error('dates separated by a ":" should have the same frequency')
            end
            ss=[ss(1:end-1),ss(end):stok];
        else
            ss=[ss,stok];
        end
    end
    x=rest;
end
end
