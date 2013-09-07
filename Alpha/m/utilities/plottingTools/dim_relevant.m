function []=dim_relevant(freq,nber_datesfunc)
if nargin<1
    freq=[];
end
if nargin<2
    [start,finish]= nber_dates(freq);
else
    nber_datesfunc=fcnchk(nber_datesfunc);
    [start,finish]= nber_datesfunc(freq);
end

% shadenber.m
%
%  Routine to shade the nberdates in a figure

curax=axis;
start_dn=[start.date_number];
finish_dn=[finish.date_number];
if start_dn(1)>curax(1)
    begin=find(start_dn>=curax(1),1,'first');  % First recession to include;
else
    begin=1;
end
if finish_dn(end)>curax(2)
    last=find(finish_dn<=curax(2),1,'first');  % First recession to include;
else
    last=numel(finish);
end

start=start(begin:last);
finish=finish(begin:last);

colorstr=[159 182 205]/256;

dim(start,finish,colorstr);

