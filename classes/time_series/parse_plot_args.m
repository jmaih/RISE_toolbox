function [this,rise_items,matlab_args]=parse_plot_args(varargin)
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

%    Defines further plotting arguments for rise time series objects 
%     The following properties allow to control for:
%     - the figure size (multiple plots): 'figsize', with default [3,3]
%     - the figure title (multiple plots): 'figtitle', with default ''
%     - the number of tick marks: 'nticks', with default 8
%     - the date format (see matlab's datestr) : 'date_format', with ''
%     - the log scale : 'logy', with default false
%     - the list of variables going to the secondary y axis : 'secondary_y'
%     - the flag for multiple plots: 'subplots' with default false
%     - vertical lines: 'vline' = '2000Q1'= '2000Q1,2003Q2' must be in the
%     same frequency as the database to be plotted
%     - horizontal lines: 'hline' =1, =[1 5.5 2]

default_nticks=plot_specs();
rise_items=struct(...
    'xrange',[],...
    'nticks',default_nticks,...
    'vline','',... vertical line
    'hline','',... horizontal line
    'logy',false,...
    'date_format',''...
    );

rise_plot_names=setdiff(fieldnames(rise_items),'xrange');
nargs=nargin;
processed=false(1,nargs);
ts_flag=false;
for iarg=1:nargs
    if processed(iarg)
        continue
    end
    viarg=varargin{iarg};
    if iarg==1 && ~isa(viarg,'ts')
        rise_items.xrange=varargin{1};
        processed(iarg)=true;
        continue
    end
    if ts_flag
        if ischar(viarg)
            viarg(isspace(viarg))=[];
            if any(strcmp(viarg,rise_plot_names))
                if iarg+1>nargs
                    error('insufficient number of arguments')
                end
                if ~isempty(varargin{iarg+1})
                    rise_items.(viarg)=varargin{iarg+1};
                end
                processed(iarg+[0,1])=true;
                continue
            end
        elseif isa(viarg,'ts')
            if numel(this)>1
                error('too many time series objects')
            end
            this=[this,{viarg}]; %#ok<AGROW>
            processed(iarg)=true;
        end
    else
        if  ~isa(viarg,'ts')
            error('time series object expected to occur in first or second position')
        end
        this={viarg};
        ts_flag=~ts_flag;
        processed(iarg)=true;
    end
end
matlab_args=varargin(~processed);

if ~isempty(rise_items.xrange)
    rise_items.xrange=process_xrange(rise_items.xrange,false);
    if numel(rise_items.xrange)<2
        error('the number of elements in the range argument cannot be less than 2')
    end
    rise_items.xrange=rise_items.xrange(1):rise_items.xrange(end);
end

if ~isempty(rise_items.hline) && ~isnumeric(rise_items.hline) && ~isscalar(rise_items.hline)
    error('argument hline should be a numeric scalar')
end

if ~isempty(rise_items.vline)
    [rise_items.vline,msg]=process_xrange(rise_items.vline,true);
    if ~isempty(msg)
        error([msg,'between vline and the database'])
    end
end

if rise_items.logy
    for its=1:numel(this)
        if any(any(any(double(this{its})<=0)))
            warning('I will have to apply the log scale despite some data being <=0')
        end
        varnames=this{its}.varnames;
        this{its}=log(this{its});
        this{its}.varnames=varnames;
    end
end

    function [xrange,msg]=process_xrange(xrange,must_be_same)
        xrange=date2serial(xrange);
        freq1=frequency2num(this{1}.frequency);
        msg='';
        if must_be_same
            [~,~,freq0]=date2year_period(xrange);
            if freq0~=freq1
                msg='frequencies are different';
            end
        else
            xrange=[xrange(1),xrange(end)];
            [year,period0,freq0]=date2year_period(xrange);
            period1=period2period(period0(:),freq0,freq1);
            xrange=dec2serial(year,period1,freq1);
        end
    end
end

    %         odd=~odd;
    %         if odd
    %             if ~ischar(viarg)
    %                 error('argument in position %0.0f is expected to be a string',iarg)
    %             end
    %             viarg(isspace(viarg))=[];
    %             if any(strcmp(viarg,rise_plot_names))
    %                 if iarg+1>nargs
    %                     error('insufficient number of arguments')
    %                 end
    %                 rise_items.(viarg)=varargin{iarg+1};
%                 processed(iarg+[0,1])=true;
%                 % double change of odd = no change
%                 continue
%             end
%         end
%         odd=~odd;