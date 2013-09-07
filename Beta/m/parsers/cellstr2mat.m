function out=cellstr2mat(CellMat) %,MatName

[eqtn_nbr,nvar]=size(CellMat);

eqtn_string='';
rows_cols=[];
for eqtn=1:eqtn_nbr
    for ivar=1:nvar
        deriv=CellMat{eqtn,ivar};
        if ~strcmp(deriv,'0')
            eqtn_string=strcat(eqtn_string,deriv,';');
            rows_cols=[rows_cols;eqtn,ivar]; %#ok<AGROW>
        end
    end
end
indices=sub2ind([eqtn_nbr,nvar],rows_cols(:,1),rows_cols(:,2));
out=struct('size',[eqtn_nbr,nvar],'string',eqtn_string(1:end-1),'indices',indices);
% brilliant: you can just take it and write it to a function

