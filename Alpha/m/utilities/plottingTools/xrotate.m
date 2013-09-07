function xrotate(angle)
if nargin==0
    angle=90;
end
tmp=gcf;
kids=get(tmp,'children');
if ischar(angle)
    angle=str2double(angle);
end
silent=true;
for ii=1:numel(kids)
    tmp=get(kids(ii),'xtick');
    testdate=serial2date(tmp(1),silent);
    if isempty(testdate)
        continue
    end
    rotateXLabels(kids(ii),angle)
end

