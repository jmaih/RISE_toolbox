function this=window(this,StartDate,EndDate,vnames,pages,sorting,trailnans)
% this=window(this,StartDate,EndDate)
if nargin<7
    trailnans=[];
    if nargin<6
        sorting=[];
        if nargin<5
            pages=[];
            if nargin<4
                vnames='';
                if nargin<3
                    EndDate=[];
                    if nargin<2
                        return
                    end
                end
            end
        end
    end
end
if isempty(sorting)
    sorting=false;
end
if isempty(trailnans)
    trailnans=true;
end
if isempty(vnames)
    vnames=this.varnames;
elseif ischar(vnames)
    vnames=cellstr(vnames);
end
if isempty(StartDate)
    StartDate=this.start;
end
if isnumeric(StartDate)
    StartDate=num2str(StartDate);
end
StartDate(isspace(StartDate))=[];
if isempty(EndDate)
    EndDate=this.finish;
end
if isnumeric(EndDate)
    EndDate=num2str(EndDate);
end
EndDate(isspace(EndDate))=[];
Data=double(this);
if isempty(pages)
    pages=1:size(Data,3);
end
pages=sort(pages);
if any(pages>size(Data,3))
    error([mfilename,':: number of pages cannot exceed ',int2str(size(Data,3))])
end
% take a robust shortcut
first_obs_date=date2serial(StartDate);
last_obs_date=date2serial(EndDate);
nobs=last_obs_date-first_obs_date+1;
Bulk=nan(nobs,this.NumberOfVariables,size(Data,3));
effective_start=max(first_obs_date,this.date_number(1));
effective_end=min(last_obs_date,this.date_number(end));
start_bulk=effective_start-first_obs_date+1;
end_bulk=effective_end-first_obs_date+1;
start_data=find(this.date_number==effective_start);
end_data=find(this.date_number==effective_end);
Bulk(start_bulk:end_bulk,:,:)=Data(start_data:end_data,:,:);
if ~isempty(vnames) && ~isempty(vnames{1})
    loc_names=locate_variables(vnames,this.varnames);
    Data=Bulk(:,loc_names,pages);
else
    Data=Bulk(:,:,pages);
end
this=rise_time_series(StartDate,Data,vnames,sorting,trailnans);
end
