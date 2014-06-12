function tank=remove_definitions(tank,definitions)
test=regexp(tank(:,1),'@{[\w/*\-+^]+}','match');%@{\w+}
test= unique([test{:}]);
out=re_evaluate(test,definitions);
for jj=1:numel(out)
    if ~isempty(out{jj})
        tank(:,1)=strrep(tank(:,1),test{jj},out{jj});
%         tank(:,1)=regexprep(tank(:,1),test{jj},out{jj});
    end
end
end
%--------------------------------------------------------------------------
function out=re_evaluate(in,definitions)
out=cell(size(in));
fields=fieldnames(definitions);
for ifield=1:numel(fields)
    ff=fields{ifield};
    tmp=definitions.(ff); %#ok<NASGU>
    eval([ff,'=tmp;'])
end
for ii=1:numel(in)
    try
        tmp=eval(in{ii}(3:end-1));
        if isa(tmp,'double')
            tmp=num2str(tmp);
        end
        out{ii}=tmp;
    end
end
end
