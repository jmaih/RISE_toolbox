classdef ts
    % ts Time series
    %
    % Constructor
    % ------------
    %
    % - [ts](ts/ts)
    %
    % Visualization
    % --------------
    %
    % - [head](ts/head)
    % - [index](ts/index)
    % - [describe](ts/describe)
    % - [display](ts/display)
    % - [jbtest](ts/jbtest)
    % - [kurtosis](ts/kurtosis)
    % - [isfinite](ts/isfinite)
    % - [isinf](ts/isinf)
    % - [isnan](ts/isnan)
    % - [ge](ts/ge)
    % - [get](ts/get)
    % - [gt](ts/gt)
    % - [le](ts/le)
    % - [lt](ts/lt)
    % - [max](ts/max)
    % - [mean](ts/mean)
    % - [median](ts/median)
    % - [min](ts/min)
    % - [mode](ts/mode)
    % - [ne](ts/ne)
    % - [quantile](ts/quantile)
    % - [range](ts/range)
    % - [skewness](ts/skewness)
    % - [sum](ts/sum)
    % - [tail](ts/tail)
    % - [var](ts/var)
    % - [std](ts/std)
    % - [spectrum](ts/spectrum)
    % - [sort](ts/sort)
    %
    % Graphing
    % ---------
    %
    % - [bar](ts/bar)
    % - [barh](ts/barh)
    % - [boxplot](ts/boxplot)
    % - [hist](ts/hist)
    % - [plot](ts/plot)
    % - [plotyy](ts/plotyy)
    %
    % Calculus
    % ---------
    %
    % - [acos](ts/acos)
    % - [acosh](ts/acosh)
    % - [acot](ts/acot)
    % - [acoth](ts/acoth)
    % - [aggregate](ts/aggregate)
    % - [allmean](ts/allmean)
    % - [apply](ts/apply)
    % - [asin](ts/asin)
    % - [asinh](ts/asinh)
    % - [atan](ts/atan)
    % - [atanh](ts/atanh)
    % - [bsxfun](ts/bsxfun)
    % - [corr](ts/corr)
    % - [corrcoef](ts/corrcoef)
    % - [cos](ts/cos)
    % - [cosh](ts/cosh)
    % - [cot](ts/cot)
    % - [coth](ts/coth)
    % - [cov](ts/cov)
    % - [cumprod](ts/cumprod)
    % - [cumsum](ts/cumsum)
    % - [decompose_series](ts/decompose_series)
    % - [eq](ts/eq)
    % - [exp](ts/exp)
    % - [hpfilter](ts/hpfilter)
    % - [interpolate](ts/interpolate)
    % - [intersect](ts/intersect)
    % - [log](ts/log)
    % - [minus](ts/minus)
    % - [mpower](ts/mpower)
    % - [mrdivide](ts/mrdivide)
    % - [mtimes](ts/mtimes)
    % - [plus](ts/plus)
    % - [power](ts/power)
    % - [rdivide](ts/rdivide)
    % - [sin](ts/sin)
    % - [sinh](ts/sinh)
    % - [transform](ts/transform)
    % - [times](ts/times)
    % - [uminus](ts/uminus)
    %
    % Lookarounds
    % ------------
    % - [pages2struct](ts/pages2struct)
    % - [subsasgn](ts/subsasgn)
    % - [subsref](ts/subsref)
    %
    % Utilities
    % ------------
    %
    % - [and](ts/and)
    % - [cat](ts/cat)
    % - [collect](ts/collect)
    % - [ctranspose](ts/ctranspose)
    % - [double](ts/double)
    % - [drop](ts/drop)
    % - [dummy](ts/dummy)
    % - [expanding](ts/expanding)
    % - [fanchart](ts/fanchart)
    % - [horzcat](ts/horzcat)
    % - [nan](ts/nan)
    % - [numel](ts/numel)
    % - [ones](ts/ones)
    % - [rand](ts/rand)
    % - [randn](ts/randn)
    % - [regress](ts/regress)
    % - [reset_start_date](ts/reset_start_date)
    % - [rolling](ts/rolling)
    % - [automatic_model_selection](ts/automatic_model_selection)
    % - [transpose](ts/transpose)
    % - [zeros](ts/zeros)
    % - [values](ts/values)
    % - [step_dummy](ts/step_dummy)
    %
    %
    % properties
    % -----------
    %
    % - [varnames] -
    % - [start] -
    % - [finish] -
    % - [frequency] -
    % - [NumberOfObservations] -
    % - [NumberOfPages] -
    % - [NumberOfVariables] -
    properties
        varnames={}
        start
        finish
        frequency
        NumberOfObservations=0;
        NumberOfPages=0;
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
        varargout=ctranspose(varargin)
        varargout=describe(varargin)
        varargout=display(varargin)
        varargout=sort(varargin)
        varargout=transpose(varargin)
        % graphing
        %---------
        varargout=bar(varargin)
        varargout=barh(varargin)
        varargout=boxplot(varargin)
        varargout=hist(varargin)
        %         varargout=line(varargin)
        varargout=plot(varargin)
        varargout=plotyy(varargin)
        % calculus
        %---------
        %         varargout=ar(varargin)
        varargout=acos(varargin)
        varargout=acosh(varargin)
        varargout=acot(varargin)
        varargout=acoth(varargin)
        varargout=allmean(varargin)
        varargout=apply(varargin)
        varargout=asin(varargin)
        varargout=asinh(varargin)
        varargout=atan(varargin)
        varargout=atanh(varargin)
        varargout=automatic_model_selection(varargin)
        varargout=bsxfun(varargin)
        varargout=corr(varargin)
        varargout=corrcoef(varargin)
        varargout=cos(varargin)
        varargout=cosh(varargin)
        varargout=cot(varargin)
        varargout=coth(varargin)
        varargout=cov(varargin)
        varargout=cumprod(varargin)
        varargout=cumsum(varargin)
        varargout=exp(varargin)
        varargout=expanding(varargin)
        varargout=fanchart(varargin)
        varargout=hpfilter(varargin)
        varargout=index(varargin)
        varargout=interpolate(varargin)
        varargout=intersect(varargin)
        varargout=jbtest(varargin)
        varargout=kurtosis(varargin)
        varargout=log(varargin)
        varargout=max(varargin)
        varargout=mean(varargin)
        varargout=median(varargin)
        varargout=min(varargin)
        varargout=minus(varargin)
        varargout=mode(varargin)
        varargout=mpower(varargin)
        varargout=mrdivide(varargin)
        varargout=mtimes(varargin)
        varargout=plus(varargin)
        varargout=power(varargin)
        varargout=rdivide(varargin)
        varargout=regress(varargin)
        varargout=rolling(varargin)
        varargout=sin(varargin)
        varargout=sinh(varargin)
        varargout=skewness(varargin)
        varargout=spectrum(varargin)
        varargout=std(varargin)
        varargout=sum(varargin)
        varargout=times(varargin)
        varargout=uminus(varargin)
        varargout=var(varargin)
        % lookarounds
        %------------
        varargout=ne(varargin)
        varargout=le(varargin)
        varargout=lt(varargin)
        varargout=gt(varargin)
        varargout=ge(varargin)
        varargout=eq(varargin)
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
        varargout=horzcat(varargin)
        varargout=isfinite(varargin)
        varargout=isinf(varargin)
        varargout=isnan(varargin)
        varargout=numel(varargin)
        varargout=pages2struct(varargin)
        varargout=quantile(varargin)
        varargout=range(varargin)
        varargout=reset_start_date(varargin)
        varargout=transform(varargin)
    end
    methods(Static)
        varargout=collect(varargin)
        varargout=dummy(varargin)
        varargout=nan(varargin)
        varargout=ones(varargin)
        varargout=rand(varargin)
        varargout=randn(varargin)
        varargout=step_dummy(varargin)
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