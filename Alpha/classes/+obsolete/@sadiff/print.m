function [c,mycall]=print(obj,wrt_index)
if nargin<2
    wrt_index=[];
end

mycall=sadiff.metastruct();
if ~isempty(wrt_index)
    indx_prefix=mycall.prefix_list{end};
end
nobj=numel(obj);
c=cell(size(obj));
wrt_i=[];
for iobj=1:nobj
    if ~isempty(wrt_index)
        wrt_i=wrt_index{iobj};
        indx=mat2str(wrt_i);
        indx(isspace(indx))=',';
        indx_name=create_handle(iobj,indx_prefix);
        [~,~,mycall]=update_line(indx,indx_name,mycall);
    end
    [c{iobj},mycall]=char(obj(iobj),wrt_i,mycall);
end
% remove unused definitions
mycall=sadiff.trim_metastruct(mycall,c);
% % format output
% fid=strcat(mycall.fid(:,1),'=',mycall.fid(:,2),';');
end
