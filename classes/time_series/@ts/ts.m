classdef ts < gogetter
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
% cumprod -   H1 line
% cumsum -   H1 line
% decompose_series -   H1 line
% describe -   H1 line
% disp -   H1 line
% double -   H1 line
% drop -   H1 line
% dummy -   H1 line
% eq -   H1 line
% exp -   H1 line
% expanding -   H1 line
% fanchart -   H1 line
% ge -   H1 line
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
% spectrum -   H1 line
% std -   H1 line
% step_dummy -   H1 line
% subsasgn -   H1 line
% subsref -   H1 line
% sum -   H1 line
% tail -   H1 line
% times -   H1 line
% transform -   H1 line
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

% frequency of the time series
frequency

end

properties(Dependent)

% number of observations in the time series
NumberOfObservations=0;

% number of pages (third dimension) of the time series
NumberOfPages=0;

% number of variables in the time series
NumberOfVariables=0;

% start time of the time series
start

% end time of the time series
finish

end

properties(Dependent,Hidden)

date_numbers

end

properties(Hidden)

data

start_date_number

description

end

methods
% constructor
%--------------
function self=ts(varargin)
% ts - constructor for time series objects
%
% ::
%
%
%   self=ts() : construct a time series with no observations
%   self=ts(start_date,data)
%   self=ts(start_date,data,varnames)
%   self=ts(start_date,data,varnames,description)
%   self=ts(start_date,data,varnames,description,trailnans)
%
% Args:
%              %
%              % - **start_date** [integer|char|serial date] : start date of
%              %   the time series. The following are admitted:
%              %   - annual data : e.g. 1990 or '1990'
%              %   - bi-annual data : e.g. '1990H1'
%              %   - Quarterly data : e.g. '1990Q3'
%              %   - monthly data : e.g. '1990M12'
%              % - **data** [numeric] : the format is nobs x nvars x npages,
%              %   where:
%              %   - **nobs** is the number of observations
%              %   - **nvars** is the number of variables
%              %   - **npages** is the number of pages (3rd dimension)
%              % - **varnames** [char|cellstr] : names of the variables in the
%              %   database
%              % - **description** [char|cellstr|{''}]: comments on each
%              %   variable in the database
%              % - **trailnans** [true|{false}]: keep or remove nans (missing
%              %   observations)
%              %
% Returns:
%    :
%              %
%              % - **self** [ts] : time series
%              %
% Note:
%              %
% Example:
            %
            % See also:
            
            self=self@gogetter();
            
            vlen=length(varargin);
            
            if vlen
                
                if isa(varargin{1},'ts')
                    
                    self=varargin{1};
                    
                else
                    
                    start_date=date2serial(varargin{1});
                    
                    [datax,vnames,dscrpt,trailnans]=dispatch_options();
                    
                    siz_data=ts.check_size(datax);
                    
                    smpl=siz_data(1);
                    
                    nvars=siz_data(2);
                    
                    do_names(vnames,'varnames')
                    
                    do_names(expanded_description(),'description',false)
                    
                    if smpl
                        
                        % throw away trailing nan observations
                        %-------------------------------------
                        first_good=1;
                        
                        last_good=smpl;
                        
                        is_valid_obs=true;
                        
                        if ~trailnans
                            
                            while all(all(isnan(datax(first_good,:,:))))
                                
                                first_good=first_good+1;
                                
                                is_valid_obs=is_valid_obs && first_good<=smpl;
                                
                                if ~is_valid_obs
                                    
                                    warning([mfilename,':: no valid observation'])
                                    
                                    break
                                    
                                end
                                
                            end
                            
                            while is_valid_obs && all(all(isnan(datax(last_good,:,:))))
                                
                                last_good=last_good-1;
                                
                            end
                            
                        end
                        
                        do_dates()
                        
                        datax=datax(first_good:last_good,:,:);
                        
                    end
                    
                    self.data=datax;
                    
                end
                
            end
            
            function out=expanded_description()
                
                out=repmat({''},1,nvars);
                
                if ischar(dscrpt)
                    
                    out(:) = {dscrpt};
                    
                else
                    
                    if numel(dscrpt)~=nvars
                        
                        if numel(dscrpt)==1
                            
                            out(:) = dscrpt;
                            
                        else
                        
                        error('number of descriptions does not match # variables')
                        
                        end
                        
                    end
                    
                    out(:) = dscrpt(:);
                    
                end
                
            end
            
            function do_names(vnames,field,check)
                
                if nargin<3,check=true; end
                
                if isempty(vnames),vnames={}; end
                
                if ischar(vnames),vnames=cellstr(vnames); end
                
                if all(cellfun(@isempty,vnames))
                    
                    vnames=repmat({''},1,nvars);
                    
                else
                    
                    if numel(vnames)~=nvars
                        
                        error(['number of columns of data should be the ',...
                            'same as the number of "',field,'"'])
                        
                    end
                    
                    if check
                        
                        for ivar=1:nvars
                            
                            thisname=vnames{ivar};
                            
                            if ~isvarname(thisname)
                                
                                warning([thisname,' is not a valid variable name'])
                                
                            end
                            
                            loc=strcmp(thisname,vnames(ivar+1:end));
                            
                            if sum(loc)>0
                                
                                warning(['variable name "',thisname,'" is duplicated'])
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
                self.(field)=vnames(:)';
                
            end
            
            function do_dates()
                
                if ~is_valid_obs
                    
                    return
                    
                end
                
                if isscalar(start_date)
                    
                    start_date=start_date+(first_good-1);
                    
                else
                    
                    start_date=start_date(first_good);
                    
                end
                
                self.start_date_number=start_date;
                
                [~,freq_]=serial2date(self.start_date_number);
                
                self.frequency=char(freq_);
                
            end
            
            function [datax,vnames,dscrpt,trailnans]=dispatch_options()
                
                datax=[];
                
                vnames='';
                
                dscrpt='';
                
                trailnans=false;
                
                if vlen>1
                    
                    datax=varargin{2};
                    
                    if vlen>2
                        
                        if ~isempty(varargin{3}),vnames=varargin{3};end
                        
                        if vlen>3
                            
                            v4=varargin{4};
                            
                            if ~isempty(v4)
                                
                                if ~(ischar(v4)||iscellstr(v4))
                                    
                                    error('dscrpt (arg # 4) must be a char or a cellstr')
                                    
                                end
                                
                                dscrpt=v4;
                                
                            end
                            
                            if vlen>4
                                
                                if ~islogical(varargin{5})
                                    
                                    error('trailnans (arg # 5) must be true or false')
                                    
                                end
                                
                                trailnans=varargin{5};
                                
                                if vlen>5
                                    
                                    error('too many input arguments')
                                    
                                end
                                
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
        function n=get.NumberOfObservations(self)
            
            n=size(self.data,1);
            
        end
        
        function n=get.NumberOfPages(self)
            
            n=size(self.data,3);
            
        end
        
        function n=get.NumberOfVariables(self)
            
            n=size(self.data,2);
            
        end
        
        function dn=get.date_numbers(self)
            
            dn=self.start_date_number;
            
            if ~isempty(dn)
                
                dn=dn+(0:self.NumberOfObservations-1);
                
            end
            
        end
        
        function s=get.start(self)
            
            sdn=self.start_date_number;
            
            if isempty(sdn)
                
                s='';
                
            else
                
                s=serial2date(sdn);
                
                s=s{1};
                
            end
            
        end
        
        function s=get.finish(self)
            
            dn=self.date_numbers;
            
            if isempty(dn)
                
                s='';
                
            else                
                
                s=serial2date(dn(end));
                
                s=s{1};
                
            end
            
        end
        
        function varargout=size(self,varargin)
            
            [varargout{1:nargout}]=size(self.data,varargin{:});
            
        end
        
        % visualization
        %--------------
        varargout=allmean(varargin)
        varargout=apply(varargin)
        varargout=bsxfun(varargin)
        varargout=chebyshev_box(varargin)
        varargout=describe(varargin)
        varargout=disp(varargin)
        varargout=expanding(varargin)
        varargout=fanchart(varargin)
        varargout=hpfilter(varargin)
        varargout=index(varargin)
        varargout=interpolate(varargin)
        varargout=intersect(varargin)
        varargout=ma_filter(varargin)
        varargout=moments(varargin)
        varargout=npdecomp(varargin)
        varargout=pdecomp(varargin)
        varargout=regress(varargin)
        varargout=rolling(varargin)
        varargout=spectrum(varargin)
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
        varargout=rmse(varargin)
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
        % utilities
        %----------
        varargout=aggregate(varargin)
        varargout=and(varargin)
        varargout=cat(varargin)
        varargout=chowlin(varargin)
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
        varargout=concatenator(varargin)
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
        varargout=check_size(varargin)
        varargout=decomp_format_output(varargin)
    end
    
    methods(Access=private)
        varargout=comparison(varargin)
        varargout=process_subs(varargin)
        varargout=ts_roll_or_expand(varargin)
    end
end