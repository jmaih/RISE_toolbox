classdef ts
    % ts Time series
    %
    % ts Methods:
    %
    % acos -   H1 line
    % acosh -   H1 line
    % acot -   H1 line
    % acoth -   H1 line
    % aggregate -   H1 line
    % allmean -   H1 line
    % and -   H1 line
    % apply -   H1 line
    % asin -   H1 line
    % asinh -   H1 line
    % atan -   H1 line
    % atanh -   H1 line
    % automatic_model_selection -   H1 line
    % bar -   H1 line
    % barh -   H1 line
    % boxplot -   H1 line
    % bsxfun -   H1 line
    % cat - concatenates time series along the specified dimension
    % collect -   H1 line
    % corr -   H1 line
    % corrcoef -   H1 line
    % cos -   H1 line
    % cosh -   H1 line
    % cot -   H1 line
    % coth -   H1 line
    % cov -   H1 line
    % ctranspose -   H1 line
    % cumprod -   H1 line
    % cumsum -   H1 line
    % decompose_series -   H1 line
    % describe -   H1 line
    % display -   H1 line
    % double -   H1 line
    % drop -   H1 line
    % dummy -   H1 line
    % eq -   H1 line
    % exp -   H1 line
    % expanding -   H1 line
    % fanchart -   H1 line
    % ge -   H1 line
    % get -   H1 line
    % gt -   H1 line
    % head -   H1 line
    % hist -   H1 line
    % horzcat -   H1 line
    % hpfilter -   H1 line
    % index -   H1 line
    % interpolate -   H1 line
    % intersect -   H1 line
    % isfinite -   H1 line
    % isinf -   H1 line
    % isnan -   H1 line
    % jbtest -   H1 line
    % kurtosis -   H1 line
    % le -   H1 line
    % log -   H1 line
    % lt -   H1 line
    % max -   H1 line
    % mean -   H1 line
    % median -   H1 line
    % min -   H1 line
    % minus -   H1 line
    % mode -   H1 line
    % mpower -   H1 line
    % mrdivide -   H1 line
    % mtimes -   H1 line
    % nan -   H1 line
    % ne -   H1 line
    % numel -   H1 line
    % ones - overloads ones for ts objects
    % pages2struct -   H1 line
    % plot -   H1 line
    % plotyy -   H1 line
    % plus -   H1 line
    % power -   H1 line
    % prctile - Percentiles of a time series (ts)
    % quantile -   H1 line
    % rand -   H1 line
    % randn -   H1 line
    % range -   H1 line
    % rdivide -   H1 line
    % regress -   H1 line
    % reset_start_date -   H1 line
    % rolling -   H1 line
    % sin -   H1 line
    % sinh -   H1 line
    % skewness -   H1 line
    % sort -   H1 line
    % spectrum -   H1 line
    % std -   H1 line
    % step_dummy -   H1 line
    % subsasgn -   H1 line
    % subsref -   H1 line
    % sum -   H1 line
    % tail -   H1 line
    % times -   H1 line
    % transform -   H1 line
    % transpose -   H1 line
    % ts - Methods:
    % uminus -   H1 line
    % values -   H1 line
    % var -   H1 line
    % zeros -   H1 line
    %
    % ts  Properties:
    %
    % varnames -   names of the variables in the database
    % start - time of the time series
    % finish -   end time of the time series
    % frequency - of the time series
    % NumberOfObservations -   number of observations in the time series
    % NumberOfPages -   number of pages (third dimension) of the time series
    % NumberOfVariables -   number of variables in the time series
    properties
        % names of the variables in the database
        varnames={}
        % start time of the time series
        start
        % end time of the time series
        finish
        % frequency of the time series
        frequency
        % number of observations in the time series
        NumberOfObservations=0;
        % number of pages (third dimension) of the time series
        NumberOfPages=0;
        % number of variables in the time series
        NumberOfVariables=0;
    end
    properties(Hidden)
        data
        date_numbers
        cell_style
    end
    methods
        % constructor
        %--------------
        function self=ts(varargin)
            % ts - constructor for time series objects
            %
            % Syntax
            % -------
            % ::
            %
            %   self=ts() : construct a time series with no observations
            %   self=ts(start_date,data)
            %   self=ts(start_date,data,varnames)
            %   self=ts(start_date,data,varnames,sorting)
            %   self=ts(start_date,data,varnames,sorting,trailnans)
            %
            % Inputs
            % -------
            %
            % - **start_date** [integer|char|serial date] : start date of
            %   the time series. The following are admitted:
            %   - annual data : e.g. 1990 or '1990'
            %   - bi-annual data : e.g. '1990H1'
            %   - Quarterly data : e.g. '1990Q3'
            %   - monthly data : e.g. '1990M12'
            % - **data** [numeric] : the format is nobs x nvars x npages,
            %   where:
            %   - **nobs** is the number of observations
            %   - **nvars** is the number of variables
            %   - **npages** is the number of pages (3rd dimension)
            % - **varnames** [char|cellstr] : names of the variables in the
            %   database
            % - **sorting** [true|{false}]: sort the columns of data
            %   according to the alphabetical order of the variable names
            % - **trailnans** [true|{false}]: keep or remove nans (missing
            %   observations)
            %
            % Outputs
            % --------
            %
            % - **self** [ts] : time series
            %
            % More About
            % ------------
            %
            % Examples
            % ---------
            %
            % See also:

            vlen=length(varargin);
            if vlen
                if isa(varargin{1},'ts')
                    self=varargin{1};
                else
                    start_date=varargin{1};
                    if ischar(start_date)
                        start_date=cellstr(start_date);
                    end
                    datax=[];if vlen>1,datax=varargin{2};end
                    vnames=''; if vlen>2,vnames=varargin{3}; end
                    sorting=[]; if vlen>3,sorting=varargin{4}; end
                    if isempty(sorting),sorting=false; end
                    trailnans=[];if vlen>4,trailnans=varargin{5}; end
                    if isempty(trailnans),trailnans=false; end
                    self.cell_style=iscell(datax);
                    if self.cell_style
                        smpl=numel(datax);
                        nvars=size(datax{1},1);
                        if nvars~=size(datax{1},2)
                            error('for matrix time series, each unit must be square')
                        end
                        npages=size(datax{1},3);
                        siz_data=[smpl,nvars,npages];
                    else
                        siz_data=size(datax);
                        smpl=siz_data(1);
                        nvars=siz_data(2);
                    end
                    if numel(siz_data)>3
                        error([mfilename,':: time series cannot have more than 3 dimensions in this environment'])
                    elseif numel(siz_data)==3
                        npages=siz_data(3);
                    else
                        npages=1;
                    end
                    if smpl
                        % throw away trailing nan observations
                        first_good=1;
                        last_good=smpl;
                        if ~trailnans && ~self.cell_style
                            while all(all(isnan(datax(first_good,:,:))))
                                first_good=first_good+1;
                                if first_good>smpl
                                    error([mfilename,':: no valid observation'])
                                end
                            end
                            while all(all(isnan(datax(last_good,:,:))))
                                last_good=last_good-1;
                            end
                        end
                        nobs=last_good-first_good+1;
                        self.NumberOfObservations=nobs;
                        self.NumberOfVariables=nvars;
                        self.NumberOfPages=npages;
                        if ~is_serial(start_date)
                            start_date=date2serial(start_date);
                        end
                        if isscalar(start_date)
                            start_date=start_date+(first_good-1)+(0:nobs-1);% <--- self.date_numbers=date2serial(start_date)+(first_good-1)+[0:nobs-1];
                        else
                            start_date=start_date(first_good:last_good);
                            if numel(start_date)~=nobs
                                error('number of input dates does not correspond to the number of observations')
                            end
                        end
                        self.date_numbers=start_date(:)';
                        [these_dates,self.frequency]=serial2date(self.date_numbers);
                        self.start=these_dates{1};
                        self.finish=these_dates{end};
                        tags=1:nvars;
                        if isempty(vnames)
                            vnames={};
                        end
                        if ischar(vnames)
                            vnames=cellstr(vnames);
                        end
                        if all(cellfun(@isempty,vnames))
                            vnames=repmat({''},1,nvars);
                        else
                            if numel(vnames)~=nvars
                                error('number of columns of data should be the same as the number of variables')
                            end
                            for ivar=1:nvars
                                thisname=vnames{ivar};
                                if ~isvarname(thisname)
                                    warning([thisname,' is not a valid variable name'])
                                end
                            end
                            if sorting
                                [vnames,tags]=sort(vnames);
                            end
                        end
                        self.varnames=vnames(:)';
                        if self.cell_style
                            if sorting
                                for iobs=1:nobs
                                    datax{iobs}=datax{iobs}(tags,tags,:); %#ok<AGROW>
                                end
                            end
                        else
                            datax=datax(first_good:last_good,tags,:);
                        end
                        self.data=datax;
                    end
                end
            end
        end
        % visualization
        %--------------
        varargout=allmean(varargin)
        varargout=apply(varargin)
        varargout=automatic_model_selection(varargin)
        varargout=bsxfun(varargin)
        varargout=chebyshev_box(varargin)
        varargout=ctranspose(varargin)
        varargout=describe(varargin)
        varargout=display(varargin)
        varargout=expanding(varargin)
        varargout=fanchart(varargin)
        varargout=hpfilter(varargin)
        varargout=index(varargin)
        varargout=interpolate(varargin)
        varargout=intersect(varargin)
        varargout=moments(varargin)
        varargout=regress(varargin)
        varargout=rolling(varargin)
        varargout=sort(varargin)
        varargout=spectrum(varargin)
        varargout=transpose(varargin)
        %         varargout=ar(varargin)
        % statistics
        %------------
        varargout=corr(varargin)
        varargout=corrcoef(varargin)
        varargout=cumprod(varargin)
        varargout=cumsum(varargin)
        varargout=jbtest(varargin)
        varargout=kurtosis(varargin)
        varargout=mean(varargin)
        varargout=median(varargin)
        varargout=mode(varargin)
        varargout=prctile(varargin)
        varargout=skewness(varargin)
        varargout=std(varargin)
        varargout=sum(varargin)
        varargout=var(varargin)
        % graphing
        %---------
        varargout=bar(varargin)
        varargout=barh(varargin)
        varargout=boxplot(varargin)
        varargout=hist(varargin)
        %         varargout=line(varargin)
        varargout=plot(varargin)
        varargout=plotyy(varargin)
        varargout=plot_real_time(varargin)
        % calculus
        %---------
        varargout=acos(varargin)
        varargout=acosh(varargin)
        varargout=acot(varargin)
        varargout=acoth(varargin)
        varargout=asin(varargin)
        varargout=asinh(varargin)
        varargout=atan(varargin)
        varargout=atanh(varargin)
        varargout=cos(varargin)
        varargout=cosh(varargin)
        varargout=cot(varargin)
        varargout=coth(varargin)
        varargout=cov(varargin)
        varargout=eq(varargin)
        varargout=exp(varargin)
        varargout=ge(varargin)
        varargout=gt(varargin)
        varargout=le(varargin)
        varargout=log(varargin)
        varargout=lt(varargin)
        varargout=max(varargin)
        varargout=min(varargin)
        varargout=minus(varargin)
        varargout=mpower(varargin)
        varargout=mrdivide(varargin)
        varargout=mtimes(varargin)
        varargout=ne(varargin)
        varargout=plus(varargin)
        varargout=power(varargin)
        varargout=rdivide(varargin)
        varargout=sin(varargin)
        varargout=sinh(varargin)
        varargout=times(varargin)
        varargout=uminus(varargin)
        % lookarounds
        %------------
        varargout=head(varargin)
        varargout=subsasgn(varargin)
        varargout=subsref(varargin)
        varargout=tail(varargin)
        varargout=values(varargin)
        varargout=double(varargin)
        varargout=get(varargin)
        % utilities
        %----------
        varargout=aggregate(varargin)
        varargout=and(varargin)
        varargout=cat(varargin)
        varargout=decompose_series(varargin)
        varargout=drop(varargin)
        varargout=dust_up(varargin)
        varargout=group(varargin)
        varargout=horzcat(varargin)
        varargout=isfinite(varargin)
        varargout=isinf(varargin)
        varargout=isnan(varargin)
        varargout=numel(varargin)
        varargout=pages2struct(varargin)
        varargout=quantile(varargin)
        varargout=range(varargin)
        varargout=reset_data(varargin)
        varargout=reset_start_date(varargin)
        varargout=transform(varargin)
    end
    methods(Static)
        varargout=collect(varargin)
        varargout=dummy(varargin)
        varargout=fold(varargin)
        varargout=nan(varargin)
        varargout=ones(varargin)
        varargout=rand(varargin)
        varargout=randn(varargin)
        varargout=step_dummy(varargin)
        varargout=unfold(varargin)
        varargout=zeros(varargin)
    end
    methods(Static,Hidden=true)
        varargout=binary_operation(varargin)
        varargout=set_locations(varargin)
        varargout=unary_operation(varargin)
    end
    methods(Access=private)
        varargout=comparison(varargin)
        varargout=main_frame(varargin)
        varargout=process_subs(varargin)
        varargout=ts_roll_or_expand(varargin)
    end
end