function this=collect(varargin)
nn=length(varargin);
exitflag= nn==0||isempty(varargin{1})||...
    (nn==1 && isa(varargin{1},'rise_time_series'));
if exitflag
    if nn==0
        this=rise_time_series.empty(0);
    else
        this=varargin{1};
    end
    return
end

cellnames=cell(1,nn);
celldata=cell(1,nn);
ii=0;
while ~isempty(varargin)
    ii=ii+1;
    db_i=varargin{1};
    if isa(db_i,'cell')
        cellnames{ii}=db_i{1};
        celldata{ii}=db_i{2};
    elseif isa(db_i,'rise_time_series')
        % if there are many names then check that nargin==1
        if db_i.NumberOfVariables>1
            if nn>1
                error([mfilename,':: in order to collect, individual databases must all have one variable'])
            end
            cellnames=db_i.varnames;
            if isempty(cellnames{1})
                error([mfilename,':: database ',int2str(ii),' must have an internal name'])
            end
            this=db_i;
            return
        else
            cellnames(ii)=db_i.varnames;
            celldata{ii}=db_i;
            if isempty(cellnames{ii})
                error([mfilename,':: database ',int2str(ii),' must have an internal name'])
            end
        end
    elseif isa(db_i,'struct')
        fields=fieldnames(db_i);
        n0=numel(fields);
        cellnames=transpose(fields);
        celldata=cell(1,n0);
        for jj=1:n0
            celldata{jj}=db_i.(cellnames{jj});
        end
        varargin=varargin(nn+1:end);
        nn=n0;
    else
        error([mfilename,':: input ',int2str(ii),' must be either instance of class ''rise_time_series'' or a cell as {''varname'',rise_time_series object} or a struct of rise_time_series objects'])
    end
    if celldata{ii}.NumberOfVariables~=1
        error([mfilename,':: database ',int2str(ii),' must have exactly one variable'])
    end
    varargin=varargin(2:end);
end
unique_names=unique(cellnames);
duplicates=false(1,numel(unique_names));
for ii=1:numel(unique_names)
    loc=find(strcmp(unique_names{ii},cellnames));
    if numel(loc)>1
        duplicates(loc)=true;
    end
end
if any(duplicates)
    disp(cellnames(duplicates))
    error([mfilename,':: the variables above are duplicated'])
end

% find the first date, the last date, the highest frequency
first_date=celldata{1}.date_number(1);
last_date=celldata{1}.date_number(end);
frequencies_codes={'','H','Q','M'
    1,2,3,4};
loc= strcmp(celldata{1}.frequency,frequencies_codes(1,:));
frequency_code=frequencies_codes{2,loc};
highest_frequency=celldata{1}.frequency;
for ii=2:nn
    if celldata{ii}.date_number(1)<first_date
        first_date=celldata{ii}.date_number(1);
    end
    if celldata{ii}.date_number(end)>last_date
        last_date=celldata{ii}.date_number(end);
    end
    loc= strcmp(celldata{ii}.frequency,frequencies_codes(1,:));
    if frequencies_codes{2,loc}>frequency_code
        highest_frequency=celldata{ii}.frequency;
        frequency_code=frequencies_codes{2,loc};
    end
end
% now that we have the start date, the end date and the
% frequency, we can go ahead and construct the new time
% series. we need to put the highest_frequency information
% into first_date and last_date somehow. The easiest way to do
% this is simply to make new dates
first_date=convert_date(first_date,highest_frequency);
last_date=convert_date(last_date,highest_frequency);

newdates=(first_date:last_date)';
newdata=nan(numel(newdates),nn);
for ii=1:nn
    % step 1, convert the dates to the new frequency
    dat_n_ii=convert_date(celldata{ii}.date_number,highest_frequency);
    oldata=double(celldata{ii});
    locations=nan(size(dat_n_ii));
    for jj=1:numel(locations)
        loc=find(newdates==dat_n_ii(jj));
        if isempty(loc)
            error([mfilename,':: date not found in new vector'])
        end
        locations(jj)=loc;
    end
    % find the date locations
    newdata(locations,ii)=oldata;
end
% trim the dates and data in case there are further trailing
% nans
first_good=find(any(~isnan(newdata),2),1,'first');
last_good=find(any(~isnan(newdata),2),1,'last');
newdates=newdates(first_good:last_good);
newdata=newdata(first_good:last_good,:);
% I could always build a special type of rise_time_series instead of
% calling rise_time_series. But the reason I need rise_time_series is because I
% started working with it earlier

% now we can safely sort the bastard
[~,tags]=sort(cellnames);
this=rise_time_series(serial2date(newdates(1)),newdata(:,tags),cellnames(tags));
end

function this=convert_date(this0,newfreq)

[~,oldfreq,year,period]=serial2date(this0);
freqtable={'','H','Q','M'
    1,2,3,4};
oldfreq_id=find(strcmp(oldfreq,freqtable(1,:)));
if isempty(oldfreq_id)
    error([mfilename,':: unknown frequency ''',oldfreq,''])
end
newfreq_id=find(strcmp(newfreq,freqtable(1,:)));
if isempty(newfreq_id)
    error([mfilename,':: unknown frequency ''',newfreq,''])
end
if oldfreq_id==newfreq_id
    this=this0;
else
    if oldfreq_id>newfreq_id
        error([mfilename,':: date conversion from high to low frequency is undefined'])
    end
    if oldfreq_id==1 % Annual
        if newfreq_id==2 % 'H'
            period=2;
        elseif newfreq_id==3 % 'Q'
            period=4;
        elseif newfreq_id==4% 'M'
            period=12;
        end
    elseif oldfreq_id==2 % H
        if newfreq_id==3% 'Q'
            period=2*period;
        elseif newfreq_id==4% 'M'
            period=6*period;
        end
    elseif oldfreq_id==3 % Q
        if newfreq_id==4% 'M'
            period=3*period;
        end
    elseif oldfreq_id==4 % M
        % this case is not expected to happen
    end
    year=num2str(year(:));
    switch newfreq
        case ''
            dat=year;
        case {'H','Q','M'}
            period=num2str(period(:));
            dat=strcat(year,newfreq,period);
    end
    this=date2serial(cellstr(dat));
end
end
