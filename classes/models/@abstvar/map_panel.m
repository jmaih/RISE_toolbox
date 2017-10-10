function [batch_rows,batch_cols]=map_panel(nvars,nx,nlags,g,ng)

% map one panel i.e. find its rows and columns in the solution

batch_rows=g:ng:ng*nvars;% batch_rows=(g-1)*nvars+1:g*nvars;

batch_cols=cell(1,nlags+1);

batch_cols{1}=g:ng:ng*nx; % batch_cols{1}=(g-1)*nx+1:g*nx;

offset=ng*nx;

for ilag=1:nlags
    
    newbatch=offset+batch_rows;
    
    batch_cols{1+ilag}=newbatch;
    
    offset=offset+nvars*ng;
    
end

batch_cols=cell2mat(batch_cols);

end