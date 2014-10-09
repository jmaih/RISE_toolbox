function [BigStart,n1_start,n2_start,n1_end,n2_end]=CombineDates(this1,this2)
if ~isequal(this1.frequency,this2.frequency)
    error([mfilename,':: databases should have the same frequency'])
end
BigStart=min(this1.date_number(1),this2.date_number(1));
BigEnd=max(this1.date_number(end),this2.date_number(end));
span=BigStart:BigEnd;
n1_start=find(this1.date_number(1)==span);
n1_end=find(this1.date_number(end)==span);
n2_start=find(this2.date_number(1)==span);
n2_end=find(this2.date_number(end)==span);
BigStart=serial2date(BigStart);
end
