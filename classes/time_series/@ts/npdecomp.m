function out=npdecomp(y0,doLog)
% NPDECOMP - Non-parametric decomposition into trend, seasonal and irregular
% components
%
% ::
%
%
%   out=NPDECOMP(y)
%   out=NPDECOMP(y,doLog)
%
% Args:
%
%    - **y** [ts] : time series to decompose
%
%    - **doLog** [true|{false}] : if log, do a multiplicative decomposition
%      otherwise the decomposition is additive
%
% Returns:
%    :
%
%    - **out** [struct] :
%      - **trend** [ts] : estimated trend
%      - **sc**    [ts] : estimated seasonal component
%      - **sa**    [ts] : seasonally adjusted data
%      - **ic**    [ts] : estimated irregular component
%
% Note:
%
%    If there are many variables and the variables are named, the first level
%    of the structure will be the names of the different variables.
%
% Example:
%
%    See also: PDECOMP
%    ---------
%

n=nargin;

check_inputs()

s=utils.time_series.freq2freq(get(y0,'frequency'));

y=double(y0);

if doLog
    
    y=log(y);
    
end

[T,nvar]=size(y);

out=struct();

% 1- first estimate of the trend component
%-----------------------------------------
trend0 = utils.time_series.mymafilter(y,6,true);

% 2- detrend the data
% -------------------
xt = y-trend0;

% 3- estimate seasonal indicator
%--------------------------------
sidx=seasonal_indices();
% S3x3 seasonal filter
sc0 = seasonal_filter(xt,3);

% 4- deseasonalize the original series
%------------------------------------
sa0 = y-sc0;

% 5- Second estimation of the trend
%------------------------------------
out.trend=henderson_filter(sa0);

% 6- detrend the data
% -------------------
xt = y-out.trend;

% 7- estimate seasonal indicator
% -------------------------------
% S3x5 seasonal filter
out.sc = seasonal_filter(xt,5);

% 8- deseasonalize the original series
% -------------------------------------
out.sa = y-out.sc;

% 9- irregular component
% -----------------------
out.ic = out.sa-out.trend;

out=ts.decomp_format_output(out,y0,doLog);

    function check_inputs()
        
        if ~isa(y0,'ts')||~get(y0,'NumberOfPages')==1
            
            error('first input must be a time series with a single page')
            
        end
        
        if n<2
            
            doLog=[];
            
        end
        
        if isempty(doLog),doLog=false;end
        
        if ~islogical(doLog)
            
            error('second input (doLog) must be a logical')
            
        end
        
    end

    function sidx=seasonal_indices()
        
        sidx=cell(1,s);
        
        for ii=1:s
            
            sidx{ii}=ii:s:T;
            
        end
        
    end

    function h13=henderson_filter(dt)
        % Henderson filter weights
        [sWH,aWH]=filter_weights('Henderson');
        
        % Apply 13-term Henderson filter
        first = 1:12;
        
        last = T-11:T;
        
        h13=dt;
        
        for jcol=1:nvar
            
            h13(:,jcol) = conv(dt(:,jcol),sWH,'same');
            
            h13(T-5:end,jcol) = conv2(dt(last,jcol),1,aWH,'valid');
            
            h13(1:6,jcol) = conv2(dt(first,jcol),1,rot90(aWH,2),'valid');
            
        end
        
    end

    function shat = seasonal_filter(xt,m)
        % S3xm seasonal filter
        [sW,aW]=filter_weights(m);
        
        % Apply filter to each period
        shat = nan(size(y));
        
        for i = 1:s
            
            nterms=m+1;
            
            ns = length(sidx{i});
            
            first = 1:nterms;
            
            last = ns-m:ns;
            
            dat = xt(sidx{i},:);
            
            for jcol=1:nvar
                
                sd = conv(dat(:,jcol),sW,'same');
                
                sd(1:nterms/2) = conv2(dat(first,jcol),1,rot90(aW,2),'valid');
                
                sd(ns-(nterms/2-1:-1:0)) = conv2(dat(last,jcol),1,aW,'valid');
                
                shat(sidx{i},jcol) = sd;
                
            end
            
        end
        
        % 13-term moving average of filtered series
        sb=utils.time_series.mymafilter(shat,6,true);
        
        shat=shat-sb;
        
    end

end

function [sW,aW]=filter_weights(m)

if ischar(m)
    
    if strcmpi(m,'henderson')
        
        sW = [-0.019;-0.028;0;.066;.147;.214;
            .24;.214;.147;.066;0;-0.028;-0.019];
        % Asymmetric weights for end of series
        aW = [
            -.034  -.017   .045   .148   .279   .421
            -.005   .051   .130   .215   .292   .353
            .061   .135   .201   .241   .254   .244
            .144   .205   .230   .216   .174   .120
            .211   .233   .208   .149   .080   .012
            .238   .210   .144   .068   .002  -.058
            .213   .146   .066   .003  -.039  -.092
            .147   .066   .004  -.025  -.042  0
            .066   .003  -.020  -.016  0      0
            .001  -.022  -.008  0      0      0
            -.026  -.011   0     0      0      0
            -.016   0      0     0      0      0    ];
        
    else
        
        error('wrong filter input')
        
    end
    
elseif m==3
    % S3x3 seasonal filter
    sW = [1/9;2/9;1/3;2/9;1/9];
    % Asymmetric weights for end of series
    aW = [.259 .407;.37 .407;.259 .185;.111 0];
    
elseif m==5
    % S3x5 seasonal filter
    % Symmetric weights
    sW = [1/15;2/15;repmat(1/5,3,1);2/15;1/15];
    % Asymmetric weights for end of series
    aW = [.150 .250 .293;
        .217 .250 .283;
        .217 .250 .283;
        .217 .183 .150;
        .133 .067    0;
        .067   0     0];
    
end

end
