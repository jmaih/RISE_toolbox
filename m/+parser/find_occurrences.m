function [occur,objectives,aliens]=find_occurrences(objectives,vlist)
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

if iscellstr(objectives)
    varlist=objectives(:)';
else
    if isa(objectives,'function_handle')
        objectives=func2str(objectives);
        if strcmp(objectives(1),'@')
            first_close=find(objectives==')',1,'first');
            objectives=objectives(first_close+1:end);
        end
    end
    if ~ischar(objectives)
        error([mfilename,':: first input must be a string or a function handle'])
    end
    if ischar(vlist)
        vlist=cellstr(vlist);
    end
    vlist=vlist(:)';
    varlist=cell2mat(strcat(vlist,'|'));
    varlist = regexp(objectives,['(?<![\w])(',varlist(1:end-1),')(?![\w])'],'match');
end
occur=ismember(vlist,varlist);

aliens=varlist(~ismember(varlist,vlist));

if ~isempty(aliens)
	warning('some variables not present in the general list')
end

end
