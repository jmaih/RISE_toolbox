function DB3=addpages(DB1,DB2)
if nargin~=2
    error([mfilename,':: number of arguments should be 2'])
end
if ~isequal(class(DB1),'rise_time_series')|| ~isequal(class(DB2),'rise_time_series')
    error([mfilename,':: both inputs must be from class rise_time_series'])
end
if isempty(DB1)
    DB3=DB2;
elseif isempty(DB2)
    DB3=DB1;
else
    if isempty(DB1.varnames{1})||isempty(DB1.varnames{2})
        error([mfilename,':: variables must have names'])
    end
    [BigStart,n1_start,n2_start,n1_end,n2_end]=CombineDates(DB1,DB2);
    if DB1.NumberOfVariables>0
        VariableNames=DB1.varnames;
        if DB2.NumberOfVariables>0
            VariableNames=[VariableNames(:);DB2.varnames(:)];
        end
    else
        VariableNames=DB2.varnames;
    end
    VariableNames=unique(VariableNames);
    nvar=numel(VariableNames);
    Data=nan(max(n1_end,n2_end),nvar,DB1.NumberOfPages+DB2.NumberOfPages);
    
    % first batch
    data1=double(DB1);
    for j=numel(VariableNames):-1:1
        vj=VariableNames{j};
        vj_id=locate_variables(vj,DB1.varnames,true);
        if ~isempty(vj_id)
            Data(n1_start:n1_end,j,1:DB1.NumberOfPages)=data1(:,vj_id,:);
        end
    end
    % second batch
    data2=double(DB2);
    for j=numel(VariableNames):-1:1
        vj=VariableNames{j};
        vj_id=locate_variables(vj,DB2.varnames,true);
        if ~isnan(vj_id)
            Data(n2_start:n2_end,j,DB1.NumberOfPages+1:end)=data2(:,vj_id,:);
        end
    end
    DB3=rise_time_series(BigStart,Data,VariableNames);
end
end
