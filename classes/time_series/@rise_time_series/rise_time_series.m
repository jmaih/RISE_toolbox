classdef rise_time_series
    properties(SetAccess = protected)
        data
        varnames
        start
        finish
        frequency
        NumberOfObservations=0;
        NumberOfPages=0;
        NumberOfVariables=0;
        date_number
    end
    methods
        varargout=addpages(varargin)
        varargout=aggregate(this,varargin)
        varargout=allmean(varargin)
        varargout=and(varargin)
        varargout=ar(varargin)
        varargout=automatic_model_selection(this,varargin)
        varargout=bar(varargin)
        varargout=bsxfun(varargin)
        varargout=corrcoef(varargin)
        varargout=cov(varargin)
        varargout=double(varargin)
        varargout=drop(varargin)
        varargout=exp(varargin)
        varargout=hist(varargin)
        varargout=horzcat(varargin)
        varargout=hpfilter(varargin)
        varargout=interpolate(varargin)
        varargout=intersect(varargin)
        varargout=jbtest(varargin)
        varargout=kurtosis(varargin)
        varargout=line(varargin)
        varargout=log(varargin)
        varargout=max(varargin)
        varargout=mode(varargin)
        varargout=mean(varargin)
        varargout=min(varargin)
        varargout=minus(varargin)
        varargout=mpower(this,varargin)
        varargout=mrdivide(varargin)
        varargout=mtimes(varargin)
        varargout=pages2struct(varargin)
        varargout=plot_separate(varargin)
        varargout=plot(varargin)
        varargout=plotyy(varargin)
        varargout=plus(varargin)
        varargout=range(varargin)
        varargout=regress(varargin)
        varargout=reset_start_date(varargin)
        varargout=rdivide(varargin)
        varargout=skewness(varargin)
        varargout=spectrum(varargin)
        varargout=std(varargin)
        varargout=subsasgn(this,varargin)
        varargout=subsref(this,varargin)
        varargout=sum(varargin)
        varargout=transform(this,varargin)
        varargout=times(varargin)
        varargout=uminus(varargin)
        varargout=var(varargin)
        varargout=window(this,varargin)
        function this=rise_time_series(StartDate,Data,varnames,sorting,trailnans)
            if nargin<5
                trailnans=false;
                if nargin<4
                    sorting=false;
                    if nargin<3
                        varnames='';
                        if nargin<2
                            Data=[];
                            if nargin<1
                                return
                            end
                        end
                    end
                end
            end
            [this.NumberOfObservations,this.NumberOfVariables,this.NumberOfPages]=size(Data);
            if this.NumberOfObservations==0;
                return
            end
            
            if ischar(StartDate)
                StartDate=cellstr(StartDate);
            end
            if numel(StartDate)>1
                if ~isequal(numel(StartDate),this.NumberOfObservations)
                    error([mfilename,':: there should be as many dates as the number of observations or just one date'])
                end
                this.date_number=date2serial(StartDate);
            else
                this.date_number=(date2serial(StartDate)+(0:this.NumberOfObservations-1))';
            end
            % throw away trailing nan observations
            smpl=size(Data,1);
            first_good=1;
            last_good=smpl;
            if ~ trailnans
                while all(all(isnan(Data(first_good,:,:))))
                    first_good=first_good+1;
                    if first_good>smpl
                        error([mfilename,':: no valid observation'])
                    end
                end
                while all(all(isnan(Data(last_good,:,:))))
                    last_good=last_good-1;
                end
            end
            this.date_number=this.date_number(first_good:last_good);
            [these_dates,this.frequency]=serial2date(this.date_number);
            Data=Data(first_good:last_good,:,:);
            % I forgot to update these
            this.NumberOfObservations=size(Data,1);
            
            this.start=these_dates{1};
            this.finish=these_dates{end};
            if ischar(varnames)
                varnames=cellstr(varnames);
            end
            [this.data,this.varnames]=CreateDataBase(these_dates,Data,varnames,sorting);
        end
        % eventually add a the following methods
        % - acf: syntax this.acf(j): autocorrelation coefficients up to j
        % - adf: syntax this.adf(lags,trendorder): augm. Dickey-Fuller order
        % - consolidate or convert a series to any frequency that would
        % result in a smaller number of observations
        % - convert: syntax this.convert(targFreq,convMethod,innerMethod,
        % ignoreMissing)
        % - moving_average
        % - summary_statistics
    end
    methods(Static)
        varargout=rand(varargin)
        varargout=collect(varargin)
        varargout=randn(varargin)
        varargout=zeros(varargin)
        varargout=nan(varargin)
        varargout=ones(varargin)
        varargout=dummy(varargin)
        varargout=step_dummy(varargin)
    end
end


function [C,varnames]=CreateDataBase(Dates,Data,varnames,sorting)
sizData=size(Data);
n=sizData(1);
k=sizData(2);
n_pages=1;
if numel(sizData)==3
    n_pages=sizData(3);
elseif numel(sizData)>3
    error([mfilename,':: time series cannot have more than 3 dimensions in this environment'])
end
if ~all(cellfun(@isempty,varnames))
    if numel(varnames)~=k
        error([mfilename,':: there should be as many variable names as the number of columns of the data'])
    end
    if sorting
        [varnames,tags]=sort(varnames);
        if ~isempty(Data)
            Data=Data(:,tags,:);
        end
    end
    Header=[{'Time'};varnames(:)];
else
    Header=[{'Time'};repmat({''},k,1)];
end
C=cell(n+1,k+1,n_pages);
for ii=1:n_pages
    C(1,:,ii)=Header';
    C(2:end,1,ii)=Dates;
    % protect against problems with num2cell when input is empty
    if ~isempty(Data)
        C(2:end,2:end,ii)=num2cell(Data(:,:,ii));
    end
end
end
