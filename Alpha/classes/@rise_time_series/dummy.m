function ts=dummy(start_date,end_date,dummy_date)

if ischar(start_date)
    start_date={start_date};
end
if ischar(dummy_date)
    dummy_date={dummy_date};
end
dummy_date=rise_date(dummy_date);
if numel(start_date)<2
    if ischar(end_date)
        end_date={end_date};
    end
    start_date=rise_date(start_date):rise_date(end_date);
end
dummy_date_numbers=[dummy_date.date_number];
sample_date_numbers=[start_date.date_number];
dummy_locs=nan(1,numel(dummy_date_numbers));
for id=1:numel(dummy_date_numbers)
    loc=find(sample_date_numbers==dummy_date_numbers(id));
    if isempty(loc)
        error(['date ',dummy_date(id).date,' could not be located in the sample'])
    end
    dummy_locs(id)=loc;
end
nobs=numel(start_date);
data=zeros(nobs,1);
data(dummy_locs)=1;

ts=rise_time_series(start_date,data);

