function list_ss=collect_symbolic_list(batch,items)
if ischar(items)
    items={items};
end
items=strcat(items,'|');
items=cell2mat(items);
items=items(1:end-1);
list_ss=regexp(batch,['(?<![\w])(',items,')[\d]+(?![\w])'],'match');
if ~isempty(list_ss) && iscell(list_ss{1})
    list_ss_={};
    for ilist=1:numel(list_ss)
        list_ss_=[list_ss_,list_ss{ilist}]; %#ok<AGROW>
    end
    list_ss=list_ss_;
end
list_ss=unique(list_ss);


% (?<!\w)(x|b|c)(\()(\d+)(\))