function [cellarray,discard_rows,Count]=trim(cellarray,min_ncount,~)
% if nargin<3
%     debug=false;
    if nargin<2
        min_ncount=3;
    end
% end
% cutoff=num2cell(1:min_ncount);

[nrows,ncols]=size(cellarray);
if ncols~=3
    error('expecting exactly 3 columns')
end
left_cells=cellarray(:,1);
right_cells=cellarray(:,2);
ncalls=cell2mat(cellarray(:,3));
trim_locations=ncalls<=min_ncount;

iter=nrows+1;
discard_rows=false(1,nrows);
% test=true;
while iter>1
    iter=iter-1;
    if ~trim_locations(iter)
        continue
    end
%     if strcmp(left_cells{iter}(1),'G')||strncmp(left_cells{iter},'if',2) % this should come from upstairs,... a restricted list or something
%         continue
%     end
%     if test
%         n_retr=ncalls(iter);
%     else
%         retrieved=regexp(right_cells(iter+1:end),['(?<!\w+)',left_cells{iter},'(?!\w+)'],'match');
%         retrieved=[retrieved{:}];
%         n_retr=numel(retrieved);
%     end
%     is_discard=true;
%     switch n_retr
%         case 0
%             if debug
%                 disp(cellarray(iter,:))
%             end
%         case cutoff
            rhs=re_parenthesize(right_cells{iter});
%             right_cells(iter+1:end)=regexprep(right_cells(iter+1:end),['(?<!\w+)',left_cells{iter},'(?!\w+)'],rhs);
            right_cells(iter+1:end)=strrep(right_cells(iter+1:end),left_cells{iter},rhs);
% %             disp(ncalls(iter))
%         otherwise
%             is_discard=false;
%     end
%     if is_discard
        left_cells(iter)=[];
        right_cells(iter)=[];
%         discard_rows(iter)=true;
%     end
end
cellarray=[left_cells,right_cells];
Count=[sum(trim_locations),nrows];
end

function rhs=re_parenthesize(rhs)
isparen=~isempty(regexp(rhs,'[\/\*\-\+]','once'));% <--[/*-+] failed on "-" and so I decided to mask them all
if isparen
    rhs=['(',rhs,')'];
end
end