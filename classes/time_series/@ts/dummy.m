function db=dummy(start_date,end_date,dummy_date)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:


start_date=date2serial(start_date);
if ischar(dummy_date)
     dummy_date=char2serial(dummy_date);
elseif isnumeric(dummy_date)
    dummy_date=date2serial(dummy_date);
end

md=max(dummy_date);
if isempty(end_date)
    start_date=start_date:md;
else
    start_date=start_date:date2serial(end_date);
end

nobs=numel(start_date);
data=zeros(nobs,1);

for id=1:numel(dummy_date)
    loc=start_date==dummy_date(id);
    if ~any(loc)
        error(['date ',serial2date(dummy_date(id)),' could not be located in the sample'])
    end
    data(loc)=1;
end

db=ts(start_date,data);
end

