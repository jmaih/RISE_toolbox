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
        TimeInfo
    end
    methods
        varargout=plot_separate(varargin)
        varargout=plot(varargin)
        varargout=plotyy(varargin)
        varargout=bar(varargin)
        varargout=line(varargin)
        varargout=mpower(this,varargin)
        function this=subsref(obj,s)
            this=obj;
            while ~isempty(s)
                if isequal(class(this),'rise_time_series')
                    if strcmp(s(1).type,'.')
                        if isequal(s(1).subs,'TimeInfo')
                            if length(s)==1 || length(s)>3
                                this=builtin('subsref',this,s);
                            else
                                if length(s)==2
                                    if isequal(s(2).subs,'date')
                                        this={this.TimeInfo.date};
                                    elseif isequal(s(2).subs,'date_number')
                                        this=[this.TimeInfo.date_number];
                                    else
                                        this=builtin('subsref',this,s);
                                    end
                                elseif length(s)==3
                                    if length(s)==3 && isequal(s(3).subs,'date_number')
                                        this=[this.TimeInfo(s(2).subs{1}).date_number];
                                    else
%                                         for kk=1:numel(s)
%                                             this=builtin('subsref',this,s(kk));
%                                         end
                                        this=builtin('subsref',this,s);
                                    end
                                end
                            end
                        else
                            this=builtin('subsref',this,s);
                        end
                        s=[];
                    else
                        access=true;
                        if isnumeric(s(1).subs{1})
                            access=false;
                            if ~isscalar(s(1).subs{1})
                                error([mfilename,':: for leads or lags, subsref must be a scalar'])
                            end
                        elseif isa(s(1).subs{1},'rise_date')
                            if ~strcmp(s(1).subs{1}(1).freq,this.frequency)
                                error([mfilename,':: date entered does not match the frequency of the reference dates'])
                            end
                            zdates={s(1).subs{1}.date};
                            locs=locate_variables(zdates,{this.TimeInfo.date},true);
                            if any(isnan(locs))
                                disp(zdates(locs))
                                error([mfilename,':: dates listed above outside the range of reference dates'])
                            end
                            by_cols=false;
                            %                 rowsid=this.TimeInfo(1).date_2_observation({s(1).subs{1}.date});
                        elseif ischar(s(1).subs{1})||iscellstr(s(1).subs{1})
                            % is it dates or is it varnames?
                            if ischar(s(1).subs{1})
                                items=cellstr(s(1).subs{1});
                            else
                                items=s(1).subs{1};
                            end
                            locs=locate_variables(items,this.varnames,true);
                            by_cols=all(~isnan(locs));
                            if ~by_cols
                                locs=locate_variables(items,{this.TimeInfo.date},true);
                                if any(isnan(locs))
                                    error([mfilename,':: variables not found in database or dates not found in date vector'])
                                end
                            end
                            %                 rowsid=this.TimeInfo(1).date_2_observation(s(1).subs{1});
                        else
                            error([mfilename,':: failed to understand this kind of subsref. Please contact junior.maih@gmail.com'])
                        end
                        datta=double(this);
                        if access
                            if by_cols
                                this=rise_time_series(this.TimeInfo,datta(:,locs,:),this.varnames(locs));
                            else
                                this=datta(locs,:,:);
                            end
                        else
                            lead_or_lag=abs(s(1).subs{1});
                            if s(1).subs{1}<0
                                datta=datta(1:end-lead_or_lag,:,:);
                                dates=this.TimeInfo(1+lead_or_lag:end);
                            elseif s(1).subs{1}>0
                                datta=datta(1+lead_or_lag:end,:,:);
                                dates=this.TimeInfo(1:end-lead_or_lag);
                            end
                            this=rise_time_series(dates,datta,this.varnames);
                        end
                        s = s(2:end); % Remove processed s from list.
                    end
                else
                    this=subsref(this,s);
                    s=[];
                end
            end
        end
        varargout=aggregate(this,varargin)
        function this=subsasgn(this,s,b)% subsasgn
            if strcmp(s(1).type,'.') || length(s)>1
                this=builtin('subsasgn',this,s,b);
            else
                vcol=1;
                rowsid=1;
                page_id=1;
                for ii=1:length(s.subs)
                    current_subs=s.subs{ii};
                    if isa(current_subs,'rise_date')
                        rowsid=locate_variables({current_subs.date},{this.TimeInfo.date},true);
                        if any(isnan(rowsid))
                            disp({current_subs(isnan(rowsid)).date})
                            error([mfilename,':: the dates above where not found in the reference date vector'])
                        end
                    elseif ischar(current_subs)
                        loc=find(strcmp(current_subs,this.varnames));
                        if ~isempty(loc)
                            vcol=loc;
                        else
                            loc=find(strcmp(current_subs,{this.TimeInfo.date}));
                            if isempty(loc)
                                error([mfilename,':: could not identify string',current_subs])
                            end
                            rowsid=loc;
                        end
                    elseif isnumeric(current_subs)
                        page_id=current_subs;
                    else
                        error([mfilename,':: I do not understand the assignment'])
                    end
                end
                datta=double(this);
                if ~isscalar(b) && (numel(b)~=numel(rowsid))
                    error([mfilename,':: number of rhs elements does not match lhs'])
                end
                datta(rowsid,vcol,page_id)=b;
                this=rise_time_series(this.TimeInfo,datta,this.varnames);
            end
        end
        function this=set_week_day(this,day)
            dn=vertcat(this.TimeInfo.date_number);
            days=rise_date.weekday(dn);
            test=unique(days);
            if numel(test)>1
                warning([mfilename,':: dates occur on different week days']) %#ok<WNTAG>
            end
            bad=days~=day;
            while any(bad)
                dn(bad)=dn(bad)+1;
                days(bad)=rise_date.weekday(dn(bad));
                bad=days~=day;
            end
            dates=rise_date.date2date(dn,this.TimeInfo(1).frequency_code);
            this=rise_time_series(dates,double(this),this.varnames);
        end
        function out=double(this)
            out=cell2mat(this.data(2:end,2:end,:));
            % cell2mat does not render the correct size when input is
            % empty. So I need to protect against that somehow as I did in
            % the past. But I will take this issue later, maybe
        end
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
            
            if this.NumberOfObservations>=1
                if isa(StartDate,'rise_date')
                    if numel(StartDate)==1
                        StartDate=StartDate+(0:this.NumberOfObservations-1);
                    end
                    if ~isequal(this.NumberOfObservations,numel(StartDate))
                        error([mfilename,':: number of observations must match the number of dates'])
                    end
                    this.TimeInfo=StartDate;
                else
                    if ischar(StartDate)
                        StartDate=cellstr(StartDate);
                    end
                    if numel(StartDate)>1
                        if ~isequal(numel(StartDate),this.NumberOfObservations)
                            error([mfilename,':: there should be as many dates as the number of observations or just one date'])
                        end
                        this.TimeInfo=rise_date(StartDate);
                    else
                        this.TimeInfo=rise_date(StartDate)+(0:this.NumberOfObservations-1);
                    end
                end
            else
                this.TimeInfo=rise_date;
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
            this.TimeInfo=this.TimeInfo(first_good:last_good);
            Data=Data(first_good:last_good,:,:);
            % I forgot to update these
            this.NumberOfObservations=size(Data,1);
            
            this.frequency=this.TimeInfo(1).freq;
            this.start=this.TimeInfo(1).date;
            this.finish=this.TimeInfo(end).date;
            if ischar(varnames)
                varnames=cellstr(varnames);
            end
            [this.data,this.varnames]=CreateDataBase({this.TimeInfo.date},Data,varnames,sorting);
        end
        function this=reset_start_date(this,startdate)
            if ~isequal(class(this),'rise_time_series')
                error([mfilename,':: input must be from class rise_time_series'])
            end
            this=rise_time_series(startdate,double(this),this.varnames);
        end
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
                trailnans=false;
            end
            if isempty(vnames)
                vnames=this.varnames;
            end
            if isempty(StartDate)
                StartDate=this.start;
            end
            if isnumeric(StartDate)
                StartDate=upper(num2str(StartDate));
            end
            if isa(StartDate,'rise_date')
                StartDate=StartDate.date;
            end
            StartDate(isspace(StartDate))=[];
            % take a robust shortcut
            startobs=find(strcmp(StartDate,{this.TimeInfo.date}));
            if isempty(startobs)
                error([mfilename,':: start date out of range of reference dates or is of a different frequency'])
            end
            if isempty(EndDate)
                EndDate=this.finish;
            end
            if isnumeric(EndDate)
                EndDate=upper(num2str(EndDate));
            end
            EndDate(isspace(EndDate))=[];
            lastobs=find(strcmp(EndDate,{this.TimeInfo.date}));
            if isempty(lastobs)
                error([mfilename,':: end date out of range of reference dates or is of a different frequency'])
            end
            loc_names=locate_variables(vnames,this.varnames);
            Data=double(this);
            if isempty(pages)
                pages=1:size(Data,3);
            end
            pages=sort(pages);
            if any(pages>size(Data,3))
                error([mfilename,':: number of pages cannot exceed ',int2str(size(Data,3))])
            end
            Data=Data(startobs:lastobs,loc_names,pages);
            this=rise_time_series(this.TimeInfo(startobs:lastobs),Data,vnames,sorting,trailnans);
        end
        function this=horzcat(varargin)
            if length(varargin)<2
                error([mfilename,':: number of arguments should be at least 2'])
            end
            for ii=1:length(varargin)
                if ~isequal(class(varargin{ii}),'rise_time_series')
                    error([mfilename,':: input ',int2str(ii),' must be from class rise_time_series'])
                end
                if ii==1
                    this=varargin{ii};
                else
                    this=this&varargin{ii};
                end
            end
        end
        function this3=and(this1,this2)
            if nargin~=2
                error([mfilename,':: number of arguments should be 2'])
            end
            if ~isequal(class(this1),'rise_time_series')|| ~isequal(class(this2),'rise_time_series')
                error([mfilename,':: both inputs must be from class rise_time_series'])
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
        function this3=plus(this1,this2)
            if isa(this1,'rise_time_series')
                if isa(this2,'rise_time_series')
                    [this1,this2]=intersect(this1,this2);
                    if ~isequal(this1.NumberOfVariables,this2.NumberOfVariables)
                        error([mfilename,':: datasets must have same number of columns'])
                    end
                    newdata=double(this1)+double(this2);
                elseif isa(this2,'double') && isscalar(this2)
                    newdata=double(this1)+this2;
                else
                    error([mfilename,':: plus operation undefined for this case'])
                end
                % Here it does not make sense to have names any more. But
                % all the same, perhaps I should have a function to rename
                % the series?
                this3=rise_time_series(this1.TimeInfo,newdata);
            elseif isa(this2,'rise_time_series') && isa(this1,'double') && isscalar(this1)
                newdata=this1+double(this2);
                % Here it does not make sense to have names any more. But
                % all the same, perhaps I should have a function to rename
                % the series?
                this3=rise_time_series(this2.TimeInfo,newdata);
            else
                error([mfilename,':: plus operation undefined for this case'])
            end
        end
        function this3=minus(this1,this2)
            if isa(this1,'rise_time_series')
                if isa(this2,'rise_time_series')
                    [this1,this2]=intersect(this1,this2);
                    if ~isequal(this1.NumberOfVariables,this2.NumberOfVariables)
                        error([mfilename,':: datasets must have same number of columns'])
                    end
                    newdata=double(this1)-double(this2);
                elseif isa(this2,'double') && isscalar(this2)
                    newdata=double(this1)-this2;
                else
                    error([mfilename,':: minus operation undefined for this case'])
                end
                this3=rise_time_series(this1.TimeInfo,newdata);
            elseif isa(this2,'rise_time_series') && isa(this1,'double') && isscalar(this1)
                newdata=this1-double(this2);
                this3=rise_time_series(this2.TimeInfo,newdata);
            else
                error([mfilename,':: minus operation undefined for this case'])
            end
        end
        function this=uminus(this)
            this=rise_time_series(this.TimeInfo,-double(this),this.varnames);
        end
        function this3=rdivide(this1,this2)
            % just to make division robust
            this3=mrdivide(this1,this2);
        end
        function this3=mrdivide(this1,this2)
            if isa(this1,'rise_time_series')
                if isa(this2,'rise_time_series')
                    [this1,this2]=intersect(this1,this2);
                    if ~isequal(this1.NumberOfVariables,this2.NumberOfVariables)
                        error([mfilename,':: datasets must have same number of columns'])
                    end
                    newdata=double(this1)./double(this2);
                elseif isa(this2,'double') && isscalar(this2)
                    newdata=double(this1)/this2;
                else
                    error([mfilename,':: divide operation undefined for this case'])
                end
                this3=rise_time_series(this1.TimeInfo,newdata);
            elseif isa(this2,'rise_time_series') && isa(this1,'double') && isscalar(this1)
                newdata=this1./double(this2);
                this3=rise_time_series(this2.TimeInfo,newdata);
            else
                error([mfilename,':: divide operation undefined for this case'])
            end
        end
        function this3=times(this1,this2)
            this3=mtimes(this1,this2);
        end
        function this3=mtimes(this1,this2)
            if isa(this1,'rise_time_series')
                if isa(this2,'rise_time_series')
                    [this1,this2]=intersect(this1,this2);
                    if ~isequal(this1.NumberOfVariables,this2.NumberOfVariables)
                        error([mfilename,':: datasets must have same number of columns'])
                    end
                    newdata=double(this1).*double(this2);
                elseif isa(this2,'double') && isscalar(this2)
                    newdata=double(this1)*this2;
                else
                    error([mfilename,':: divide operation undefined for this case'])
                end
                this3=rise_time_series(this1.TimeInfo,newdata);
            elseif isa(this2,'rise_time_series') && isa(this1,'double') && isscalar(this1)
                newdata=this1*double(this2);
                this3=rise_time_series(this2.TimeInfo,newdata);
            else
                error([mfilename,':: divide operation undefined for this case'])
            end
        end
        function this=log(this,keep_name)
            if nargin<2
                keep_name=false;
            end
            if keep_name
                this=rise_time_series(this.TimeInfo,log(double(this)),this.varnames);
            else
                this=rise_time_series(this.TimeInfo,log(double(this)));
            end
        end
        function this=exp(this)
            % Here it does not make sense to have names any more. But
            % all the same, perhaps I should have a function to rename
            % the series?
            this=rise_time_series(this.TimeInfo,exp(double(this)));
        end
        function [this1,this2]=intersect(this1,this2)
            if nargin~=2
                error([mfilename,':: function must have 2 arguments'])
            end
            if ~isa(this1,'rise_time_series')||~isa(this2,'rise_time_series')
                error([mfilename,':: both arguments should be time series'])
            end
            if ~strcmp(this1.TimeInfo(1).freq,this2.TimeInfo(1).freq)
                error([mfilename,':: datasets must have same frequency'])
            end
            [~,I1,I2] = intersect([this1.TimeInfo.date_number],[this2.TimeInfo.date_number]);
            if isempty(I1)||isempty(I2)
                error([mfilename,':: don''t have common dates'])
            end
            C=this1.TimeInfo(I1);
            vname1=this1.varnames;
            vname2=this2.varnames;
            this1=double(this1);
            this2=double(this2);
            this1=rise_time_series(C,this1(I1,:),vname1);
            this2=rise_time_series(C,this2(I2,:),vname2);
        end
        function m=mean(this)
            this=double(this);
            n=size(this,2);
            if size(this,3)>1
                error([mfilename,':: this operation is only defined for databases with one page'])
            end
            m=nan(1,n);
            for ii=1:n
                dd=this(:,ii);
                dd=dd(~isnan(dd));
                if isempty(dd)
                    error([mfilename,':: no valid observations to compute the mean in column ',int2str(ii)])
                end
                m(ii)=mean(dd);
            end
        end
        function m=min(this)
            this=double(this);
            n=size(this,2);
            if size(this,3)>1
                error([mfilename,':: this operation is only defined for databases with one page'])
            end
            m=nan(1,n);
            for ii=1:n
                dd=this(:,ii);
                dd=dd(~isnan(dd));
                if isempty(dd)
                    error([mfilename,':: no valid observations to compute the min in column ',int2str(ii)])
                end
                m(ii)=min(dd);
            end
        end
        function m=max(this)
            this=double(this);
            n=size(this,2);
            if size(this,3)>1
                error([mfilename,':: this operation is only defined for databases with one page'])
            end
            m=nan(1,n);
            for ii=1:n
                dd=this(:,ii);
                dd=dd(~isnan(dd));
                if isempty(dd)
                    error([mfilename,':: no valid observations to compute the max in column ',int2str(ii)])
                end
                m(ii)=max(dd);
            end
        end
        function m=sum(this)
            this=double(this);
            n=size(this,2);
            if size(this,3)>1
                error([mfilename,':: this operation is only defined for databases with one page'])
            end
            m=nan(1,n);
            for ii=1:n
                dd=this(:,ii);
                dd=dd(~isnan(dd));
                if isempty(dd)
                    error([mfilename,':: no valid observations to compute the sum in column ',int2str(ii)])
                end
                m(ii)=sum(dd);
            end
        end
        function m=range(this)
            this=double(this);
            n=size(this,2);
            if size(this,3)>1
                error([mfilename,':: this operation is only defined for databases with one page'])
            end
            m=nan(2,n);
            for ii=1:n
                dd=this(:,ii);
                dd=dd(~isnan(dd));
                if isempty(dd)
                    error([mfilename,':: no valid observations to compute the range in column ',int2str(ii)])
                end
                m(1,ii)=min(dd);
                m(2,ii)=max(dd);
            end
        end
        function m=std(this,varargin)
            this=double(this);
            n=size(this,2);
            if size(this,3)>1
                error([mfilename,':: this operation is only defined for databases with one page'])
            end
            m=nan(1,n);
            for ii=1:n
                dd=this(:,ii);
                dd=dd(~isnan(dd));
                if isempty(dd)
                    error([mfilename,':: no valid observations to compute the standard deviation in column ',int2str(ii)])
                end
                m(ii)=std(dd,varargin{:});
            end
        end
        function m=cov(this,varargin)
            if ~isempty(varargin) && isa(varargin{1},'rise_time_series')
                this=this & varargin{1};
                varargin=varargin(2:end);
            end
            this=double(this);
            if size(this,3)>1
                error([mfilename,':: this operation is only defined for databases with one page'])
            end
            nanrows=any(isnan(this),2);
            dd=this(~nanrows,:);
            if isempty(dd)
                error([mfilename,':: no valid observations to compute the covariance '])
            end
            m=cov(dd,varargin{:});
        end
        function m=var(this,varargin)
            this=double(this);
            n=size(this,2);
            if size(this,3)>1
                error([mfilename,':: this operation is only defined for databases with one page'])
            end
            m=nan(1,n);
            for ii=1:n
                dd=this(:,ii);
                dd=dd(~isnan(dd));
                if isempty(dd)
                    error([mfilename,':: no valid observations to compute the standard deviation in column ',int2str(ii)])
                end
                m(ii)=var(dd,varargin{:});
            end
        end
        function varargout=corrcoef(this,varargin)
            if ~isempty(varargin) && isa(varargin{1},'rise_time_series')
                this=this & varargin{1};
                varargin=varargin(2:end);
            end
            this=double(this);
            if size(this,3)>1
                error([mfilename,':: this operation is only defined for databases with one page'])
            end
            nanrows=any(isnan(this),2);
            dd=this(~nanrows,:);
            if isempty(dd)
                error([mfilename,':: no valid observations to compute the correlation '])
            end
            [R,P,RLO,RUP]=corrcoef(dd,varargin{:});
            nout=nargout;
            varargout{1}=R;
            if nout>1
                varargout{2}=P;
                if nout>2
                    varargout{3}=RLO;
                    if nout>3
                        varargout{4}=RUP;
                    end
                end
            end
        end
        function K=kurtosis(this,varargin)
            this=double(this);
            if size(this,3)>1
                error([mfilename,':: this operation is only defined for databases with one page'])
            end
            K=kurtosis(this,varargin{:});
        end
        function S=skewness(this,varargin)
            this=double(this);
            if size(this,3)>1
                error([mfilename,':: this operation is only defined for databases with one page'])
            end
            S=skewness(this,varargin{:});
        end
        function varagout=mode(this,varargin)
            this=double(this);
            if size(this,3)>1
                error([mfilename,':: this operation is only defined for databases with one page'])
            end
            nanrows=any(isnan(this),2);
            dd=this(~nanrows,:);
            if isempty(dd)
                error([mfilename,':: no valid observations to compute the mode '])
            end
            varagout=mode(dd,varargin{:});
        end
        function varagout=hist(this,varargin)
            this=double(this);
            if size(this,3)>1
                error([mfilename,':: this operation is only defined for databases with one page'])
            end
            nanrows=any(isnan(this),2);
            dd=this(~nanrows,:);
            if isempty(dd)
                error([mfilename,':: no valid observations to compute the mode '])
            end
            varagout=hist(dd,varargin{:});
        end
        function m=allmean(this)
            vnames=this.varnames;
            this=double(this);
            if size(this,3)>1
                error([mfilename,':: this operation is only defined for databases with one page'])
            end
            n=size(this,2);
            m=cell(5,n+1);
            m(2:end,1)={'Harmonic','Geometric','Arithmetic','Quadratic'};
            if all(~cellfun(@isempty,vnames))
                m(1,2:end)=vnames;
            end
            r=[-1,0,1,2]+eps;
            for ii=1:n
                dd=this(:,ii);
                dd=dd(~isnan(dd));
                if isempty(dd)
                    error([mfilename,':: no valid observations to compute the mean in column ',int2str(ii)])
                end
                nd=numel(dd);
                for jj=1:4
                    m{jj+1,ii+1}=(sum(dd.^r(jj))/nd)^(1/r(jj));
                end
            end
        end
        function [h,p,jbstat,critval]=jbtest(this,varargin)
            % -H: H=0 means null hypothesis ("the data are normally
            % distributed") cannot be rejected at the 5% significance
            % level. H=1 means null can be rejected at the 5% level.
            % - P: p-value
            % - JBSTAT: test statistic
            % - CRITVAL: critical value for test
            % JBTEST treats NaNs in X as missing values, and ignores them.
            % Test Statistic:  JBSTAT = N*(SKEWNESS^2/6 + (KURTOSIS-3)^2/24),
            %                   where N is the sample size and the kurtosis of
            %                   the normal distribution is defined as 3.
            %  H = JBTEST(X,ALPHA) performs the test at significance level
            %  ALPHA
            n=numel(this);
            this=double(this);
            if size(this,3)>1
                error([mfilename,':: this operation is only defined for databases with one page'])
            end
            h=nan(1,n);
            p=nan(1,n);
            jbstat=nan(1,n);
            critval=nan(1,n);
            for ii=1:n
                dd=this(:,ii);
                dd=dd(~isnan(dd));
                if isempty(dd)
                    error([mfilename,':: no valid observations to compute the standard deviation in column ',int2str(ii)])
                end
                [h(ii),p(ii),jbstat(ii),critval(ii)]=jbtest(dd,varargin{:});
            end
        end
        function [sw,jj,T]=spectrum(this)
            y=double(this);
            if size(y,3)>1
                error([mfilename,':: this operation is only defined for databases with one page'])
            end
            bad=any(isnan(y),2);
            y=y(~bad,:);
            if isempty(y)
                error([mfilename,':: no valid data to construct the periodogram'])
            end
            [T,n]=size(y);
            
            % autocovariances
            y=bsxfun(@minus,y,mean(y,1));
            gam=nan(T,n);
            for ii=1:T
                gam(ii,:)=sum(y(ii:T,:).*y(1:T-ii+1,:),1);
            end
            gam=gam/T;
            % spectrum
            jj=(1:T-1)';
            w=2*pi*jj/T;
            sw=nan(T-1,n);
            for ii=1:numel(w)
                cwj=cos(w(ii)*jj);
                sw(ii,:)=1/(2*pi)*(gam(1,:)+2*sum(bsxfun(@times,gam(2:end,:),cwj),1));
            end
        end
        function this=ones(this)
            % Adds a column of ones at the end of the series
            datta=double(this);
            if size(datta,3)>1
                error([mfilename,':: this operation is only defined for databases with one page'])
            end
            smpl=size(datta,1);
            datta=[datta,ones(smpl,1)];
            this=rise_time_series(this.TimeInfo,datta);
        end
        function [B,BINT,R,RINT,STATS]=regress(this,this2,varargin)
            % - B: vector of regression coefficients in the linear model Y = X*B.
            % - BINT: of 95% confidence intervals for B
            % - R: vector of residuals
            % - RINT: matrix of intervals that can be used to diagnose
            % outliers.  If RINT(i,:) does not contain zero, then the i-th
            % residual is larger than would be expected, at the 5%
            % significance level.  This is evidence that the I-th
            % observation is an  outlier.
            % - STATS: vector containing, in the following order, the R-square
            % statistic, the F statistic and p value for the full model,
            % and an estimate of the error variance.
            % example: y=rise_time_series(1990,rand(100,1)); % random series
            % X=y(-1)&y(-2)&y(-3); % columns of lags
            % X=ones(X); % add a column of ones
            % [B,BINT,R,RINT,STATS]=regress(y,X)
            if isa(this,'rise_time_series') && isa(this2,'rise_time_series')
                [this,this2]=intersect(this,this2);
                y=double(this);
                if size(y,3)>1
                    error([mfilename,':: this operation is only defined for databases with one page'])
                end
                if size(y,2)>1
                    error([mfilename,':: first argument must have only one column of data'])
                end
                X=double(this2);
                if size(X,3)>1
                    error([mfilename,':: this operation is only defined for databases with one page'])
                end
                [B,BINT,R,RINT,STATS]=regress(y,X,varargin{:});
            else
                error([mfilename,':: both arguments should be ''rise_time_series'' objects'])
            end
        end
        function [B,BINT,R,RINT,STATS]=ar(this,lag,const)
            % - B: vector of regression coefficients in the linear model Y = X*B.
            % - BINT: of 95% confidence intervals for B
            % - R: vector of residuals
            % - RINT: matrix of intervals that can be used to diagnose
            % outliers.  If RINT(i,:) does not contain zero, then the i-th
            % residual is larger than would be expected, at the 5%
            % significance level.  This is evidence that the I-th
            % observation is an  outlier.
            % - STATS: vector containing, in the following order, the R-square
            % statistic, the F statistic and p value for the full model,
            % and an estimate of the error variance.
            if nargin<3
                const=true;
                if nargin<2
                    lag=1;
                end
            end
            Z=double(this);
            if size(Z,3)>1
                error([mfilename,':: this operation is only defined for databases with one page'])
            end
            if size(Z,2)>1
                error([mfilename,':: first argument must have only one column of data'])
            end
            y=Z(lag+1:end);
            smpl=size(y,1);
            X=nan(smpl,lag);
            for ii=1:lag
                X(:,ii)=Z((lag+1:end)-ii);
            end
            if const
                X=[X,ones(smpl,1)];
            end
            [B,BINT,R,RINT,STATS]=regress(y,X);
        end
        function this=pages2struct(this0)
            vnames=this0.varnames;
            if numel(unique(vnames))~=this0.NumberOfVariables
                error([mfilename,':: number of unique variable names different from number of columns of data matrix'])
            end
            this=struct();
            datta=double(this0);
            for ii=1:this0.NumberOfVariables
                this.(this0.varnames{ii})=rise_time_series(this0.TimeInfo,permute(datta(:,ii,:),[1,3,2]));
            end
        end
        function this=interpolate(this,method)
            if nargin<2
                method='spline';
            end
            tmp=double(this);
            [n,p]=size(tmp);
            if p>1
                error([mfilename,':: routine written to handle time series with only one variable.'])
            end
            w=(1:n)';
            nanrows=isnan(tmp);
            ti=w(~nanrows);
            yi=tmp(~nanrows);
            nx=numel(ti);
            wnan=w(nanrows);
            switch lower(method)
                case 'lagrange'
                    fnan=0;
                    for k=1:nx
                        p=1;
                        for j=[1:k-1,k+1:nx]
                            p=p.*(wnan-ti(j))/(ti(k)-ti(j));
                        end
                        fnan=fnan+p*yi(k);
                    end
                case 'spline'
                    h=nan(nx-1,1);
                    b=nan(nx-1,1);
                    u=nan(nx-1,1);
                    v=nan(nx-1,1);
                    z=nan(nx,1);
                    for ii=1:nx-1
                        h(ii)=ti(ii+1)-ti(ii);
                        b(ii)=6*(yi(ii+1)-yi(ii))/h(ii);
                    end
                    u(1)=2*(h(1)+h(2));
                    v(1)=b(2)-b(1);
                    for ii=2:nx-1
                        u(ii)=2*(h(ii)+h(ii-1))-h(ii-1)^2/u(ii-1);
                        v(ii)=b(ii)-b(ii-1)-h(ii-1)*v(ii-1)/u(ii-1);
                    end
                    z(nx)=0;
                    for ii=nx-1:-1:2
                        z(ii)=(v(ii)-h(ii)*z(ii+1))/u(ii);
                    end
                    z(1)=0;
                    fnan=nan(size(wnan));
                    for ii=1:nx-1
                        inan=wnan>ti(ii) & wnan<ti(ii+1);
                        fnan(inan)=z(ii)/(6*h(ii))*(ti(ii+1)-wnan(inan)).^3+...
                            z(ii+1)/(6*h(ii))*(wnan(inan)-ti(ii)).^3+...
                            (yi(ii+1)/h(ii)-z(ii+1)*h(ii)/6)*(wnan(inan)-ti(ii))+...
                            (yi(ii)/h(ii)-z(ii)*h(ii)/6)*(ti(ii+1)-wnan(inan));
                    end
                    tmp(nanrows)=fnan;
                    this=rise_time_series(this.TimeInfo,tmp);
                otherwise
                    error([mfilename,':: other interpolation methods not implemented yet'])
            end
        end
        function DB3=addpages(DB1,DB2)
            if nargin~=2
                error([mfilename,':: number of arguments should be 2'])
            end
            if ~isequal(class(DB1),'rise_time_series')|| ~isequal(class(DB2),'rise_time_series')
                error([mfilename,':: both inputs must be from class rise_time_series'])
            end
            if isempty(DB1)
                DB3=DB2;
            elseif isempty(DB2)
                DB3=DB1;
            else
                if isempty(DB1.varnames{1})||isempty(DB1.varnames{2})
                    error([mfilename,':: variables must have names'])
                end
                [BigStart,n1_start,n2_start,n1_end,n2_end]=CombineDates(DB1,DB2);
                if DB1.NumberOfVariables>0
                    VariableNames=DB1.varnames;
                    if DB2.NumberOfVariables>0
                        VariableNames=[VariableNames(:);DB2.varnames(:)];
                    end
                else
                    VariableNames=DB2.varnames;
                end
                VariableNames=unique(VariableNames);
                nvar=numel(VariableNames);
                Data=nan(max(n1_end,n2_end),nvar,DB1.NumberOfPages+DB2.NumberOfPages);
                
                % first batch
                data1=double(DB1);
                for j=numel(VariableNames):-1:1
                    vj=VariableNames{j};
                    vj_id=locate_variables(vj,DB1.varnames,true);
                    if ~isempty(vj_id)
                        Data(n1_start:n1_end,j,1:DB1.NumberOfPages)=data1(:,vj_id,:);
                    end
                end
                % second batch
                data2=double(DB2);
                for j=numel(VariableNames):-1:1
                    vj=VariableNames{j};
                    vj_id=locate_variables(vj,DB2.varnames,true);
                    if ~isnan(vj_id)
                        Data(n2_start:n2_end,j,DB1.NumberOfPages+1:end)=data2(:,vj_id,:);
                    end
                end
                DB3=rise_time_series(BigStart,Data,VariableNames);
            end
        end
        function this=drop(this,varargin)
            survive=true(1,this.NumberOfVariables);
            for ii=1:length(varargin)
                ids=locate_variables(varargin{ii},this.varnames);
                survive(ids)=false;
            end
            tmp_data=double(this);
            this=rise_time_series(this.TimeInfo,tmp_data(:,survive),this.varnames(survive));
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
        % - zeros
        varargout=automatic_model_selection(this,varargin)
    end
    methods(Static)
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
            first_date=celldata{1}.TimeInfo(1);
            last_date=celldata{1}.TimeInfo(end);
            highest_frequency=celldata{1}.TimeInfo(1).frequency_code;
            new_wknd_days=celldata{1}.TimeInfo(1).week_end_days;
            new_weekly_day=rise_date.weekday(celldata{1}.TimeInfo(1).date_number);
            for ii=2:nn
                if celldata{ii}.TimeInfo(1)<first_date
                    first_date=celldata{ii}.TimeInfo(1);
                end
                if celldata{ii}.TimeInfo(end)>last_date
                    last_date=celldata{ii}.TimeInfo(end);
                end
                if celldata{ii}.TimeInfo(1).frequency_code>highest_frequency
                    highest_frequency=celldata{ii}.TimeInfo(1).frequency_code;
                    new_wknd_days=celldata{ii}.TimeInfo(1).week_end_days;
                    new_weekly_day=rise_date.weekday(celldata{ii}.TimeInfo(1).last_business_day);
                end
            end
            % now that we have the start date, the end date and the
            % frequency, we can go ahead and construct the new time
            % series. we need to put the highest_frequency information
            % into first_date and last_date somehow. The easiest way to do
            % this is simply to make new dates
            first_date=first_date.convert_date(highest_frequency,new_wknd_days,new_weekly_day);
            last_date=last_date.convert_date(highest_frequency,new_wknd_days,new_weekly_day);
            
            newdates=first_date:last_date;
            
            newdata=nan(numel(newdates),nn);
            new_date_numbers=[newdates.date_number];
            for ii=1:nn
                % step 1, convert the dates to the new frequency
                dat_n_ii=convert_date(celldata{ii}.TimeInfo,highest_frequency,new_wknd_days,new_weekly_day);
                dat_n_ii=[dat_n_ii.date_number];
                oldata=double(celldata{ii});
                locations=nan(size(dat_n_ii));
                for jj=1:numel(locations)
                    loc=find(new_date_numbers==dat_n_ii(jj));
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
            this=rise_time_series(newdates,newdata(:,tags),cellnames(tags));
        end
        %         function SaveDataBase(this,SaveUnderName,extension)
        %             % SaveDataBase(this,SaveUnderName,extension)
        %             if nargin<3
        %                 extension='mat';
        %                 if nargin<2
        %                     error([mfilename,':: at least the database and the name to save under should be provided'])
        %                 end
        %             end
        %
        %             if ~isa(this,'rise_time_series')
        %                 error([mfilename,':: first input must be an object of class rise_time_series'])
        %             end
        %
        %             switch extension
        %                 case 'xls'
        %                     xlswrite(SaveUnderName,this.data(:,:,1))
        %                 case 'mat'
        %                     datatosave=cell2mat(this.data(2:end,2:end,1));
        %                     for j=this.NumberOfVariables:-1:1
        %                         vj=deblank(this.varnames(j,:));
        %                         DataStructure.(vj)=datatosave(:,j); %#ok<STRNU>
        %                     end
        %                     eval(['save ',SaveUnderName,' -struct DataStructure'])
        %                 case 'm'
        %                     datatosave=cell2mat(this.data(2:end,2:end,1));
        %                     fid=fopen([SaveUnderName,'.m'],'W');
        %                     fprintf(fid,'%s \n\n',['% dataset automatically generated on ',...
        %                         datestr(now)]);
        %                     fprintf(fid,'%s \n\n',['% running from ',this.start,' to ',...
        %                         this.finish,' (',int2str(this.NumberOfObservations),' observations)  ']);
        %                     fprintf(fid,'%s \n\n',['% ',int2str(this.NumberOfVariables),...
        %                         ' variables. More details below']);
        %                     fprintf(fid,'%s \n','datamat=[');
        %                     format=' \n';
        %                     for j=this.NumberOfVariables:-1:1
        %                         format=[' %12.8f',format]; %#ok<AGROW>
        %                     end
        %                     fprintf(fid,format,transpose(datatosave));
        %                     fprintf(fid,'%s \n','];');
        %                     for j=this.NumberOfVariables:-1:1
        %                         fprintf(fid,'%s \n',[this.varnames(j,:),'=datamat(:,',int2str(j),');']);
        %                     end
        %                     fclose(fid);
        %                 otherwise
        %                     error([mfilename,':: format ',extension,' unknown'])
        %             end
        %         end
    end
end

function [BigStart,n1_start,n2_start,n1_end,n2_end]=CombineDates(this1,this2)
BigStart=min(this1.TimeInfo(1),this2.TimeInfo(1));
if ~isequal(this1.TimeInfo(1).freq,this1.TimeInfo(2).freq)
    error([mfilename,':: databases should have the same frequency'])
end
n1_start=BigStart.date_2_observation(this1.TimeInfo(1).date);
n1_end=BigStart.date_2_observation(this1.TimeInfo(end).date);
n2_start=BigStart.date_2_observation(this2.TimeInfo(1).date);
n2_end=BigStart.date_2_observation(this2.TimeInfo(end).date);
BigStart=BigStart.date;
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
