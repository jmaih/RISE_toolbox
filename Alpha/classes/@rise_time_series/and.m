function this3=and(this1,this2)
if nargin~=2
    error([mfilename,':: number of arguments should be 2'])
end
if isempty(this2)
    this3=this1;
    return
end
if isempty(this1)
    this3=this2;
    return
end
[BigStart,n1_start,n2_start,n1_end,n2_end]=CombineDates(this1,this2);
nvar=this1.NumberOfVariables+this2.NumberOfVariables;
Datta1=double(this1);
Datta2=double(this2);

Data=nan(max(n1_end,n2_end),nvar);
if this1.NumberOfVariables
    Data(n1_start:n1_end,1:this1.NumberOfVariables)=Datta1;
end
Data(n2_start:n2_end,this1.NumberOfVariables+1:end)=Datta2;
this3=rise_time_series(BigStart,Data,[this1.varnames(:);this2.varnames(:)]);
end