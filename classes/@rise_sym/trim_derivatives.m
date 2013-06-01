function [cellarray,discard_rows,Count]=trim_derivatives(cellarray,debug)
if nargin<2
    debug=false;
end

[nrows,ncols]=size(cellarray);
if ncols~=2
    error('expecting exactly 2 columns')
end
left_cells=cellarray(:,1);
right_cells=cellarray(:,2);
iter=nrows+1;
discard_rows=false(1,nrows);
while iter>1
    iter=iter-1;
    if strcmp(left_cells{iter}(1),'G')
        continue
    end
    retrieved=regexp(right_cells(iter+1:end),['(?!w)',left_cells{iter},'(?<!w)'],'match');
    retrieved=[retrieved{:}];
    if isempty(retrieved)
        if debug
            disp(cellarray(iter,:))
        end
        left_cells(iter)=[];
        right_cells(iter)=[];
        discard_rows(iter)=true;
    end
end
cellarray=[left_cells,right_cells];
Count=[sum(discard_rows),nrows];