function [data,varinfo] = haver2rise(tickers)
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

%% haver2rise.m

yourPathToHaver = '//imfdata/econ/DATA/DLX/DATA/';

% tickers={Names,PSEUDO,FREQ}
nrows=size(tickers,1);

haver_series_names=cell(nrows,1);
haver_series_databases=cell(nrows,1);
haver_series_periods=cell(nrows,1);
model_names=tickers(:,2);
varinfo=cell(2,nrows);
varinfo(1,:)=model_names;

for ii=1:nrows
    vi=tickers{ii,1};
    ampers=strfind(vi,'@');
    if isempty(ampers)
        error([mfilename,':: wrong ticker name ',vi])
    end
    haver_series_names{ii}=vi(1:ampers-1);
    haver_series_databases{ii}=vi(ampers+1:end);
    if size(tickers,2)>2
        haver_series_periods{ii}=upper(tickers{ii,3});
    else
        haver_series_periods{ii}=nan;
    end
end
databases=unique(haver_series_databases);
for db=1:numel(databases)
    locs=strcmp(databases{db},haver_series_databases);
    summary.(databases{db})=haver_series_names(locs);
    summary_periods.(databases{db})=haver_series_periods(locs);
end

%% Get data

data=struct();
for ii = 1 : numel(databases)
    try
        fprintf('\n Attempting to access %s ... \n',databases{ii})
        connect = haver([yourPathToHaver,databases{ii},'.dat']); %create connection object
        fprintf('Connection to the HAVER database %s established. \n',databases{ii})
    catch me1
        fprintf('MATLAB SAYS: \n %s \n',me1.message) %if you catch an exception, then display it ...
        continue % and continue to next iteration
    end
    varnames=summary.(databases{ii});
    Frequence=summary_periods.(databases{ii});
    
    for v=1:numel(varnames)
        success=false;
        try
            this_info=info(connect,varnames(v));
            this_data=fetch(connect,varnames(v)); % [date_numbers,values]
            success=true;
        end
        if ~success
            warning([mfilename,':: unable to find variable ',varnames{v}])
        else
            switch this_info.Frequency
                case 'A'
                    this_dates=datestr(this_data(:,1),'YYYY');
                    this_dates=rise_date(this_dates);
                case 'Q'
                    this_dates=datestr(this_data(:,1),'QQ-YYYY');
                    years=this_dates(:,4:end);
                    periods=this_dates(:,2);
                    this_dates=rise_date(strcat(years,'Q',periods));
                case 'M'
                    this_dates=datestr(this_data(:,1),'mm/dd/yyyy');
                    years=this_dates(:,7:end);
                    periods=this_dates(:,1:2);
                    this_dates=rise_date(strcat(years,'M',periods));
                otherwise
                    error([mfilename,':: frequency ',this_info.Frequency,' not implemented yet'])
            end
            vloc=strcmp(varnames(v),haver_series_names);
            data.(model_names{vloc})=ts(this_dates,this_data(:,2),model_names{vloc});
            oldFreq=data.(model_names{vloc}).frequency;
            if ~isnan(Frequence{v}) && ~isequal(Frequence{v},oldFreq)
                data.(model_names{vloc})=aggregate(data.(model_names{vloc}),Frequence{v});
                warning([mfilename,':: series ',model_names{vloc},' has been aggregated from ',oldFreq,' to ',Frequence{v},' frequency'])
            end
            varinfo{2,vloc}=this_info;
        end
    end
end
