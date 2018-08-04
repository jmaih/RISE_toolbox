function out=pdecomp(y0,doLog,dorder)
% PDECOMP - Parametric decomposition into trend, seasonal and irregular
% components
%
% ::
%
%
%   out=PDECOMP(y)
%   out=PDECOMP(y,doLog)
%   out=PDECOMP(y,doLog,dorder)
%
% Args:
%
%    - **y** [ts] : time series to decompose
%
%    - **doLog** [true|{false}] : if log, do a multiplicative decomposition
%      otherwise the decomposition is additive
%
%    - **dorder** [integer|{2}] : detrending order
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
%    See also: NPDECOMP
%    ---------
%

n=nargin;

check_inputs()

y=double(y0);

if doLog
    
    y=log(y);
    
end

[T,nvar]=size(y);

out=struct();

% parametric trend
%-----------------
out.trend=parametric_trend();

% detrend the data
% ----------------
xt=y-out.trend;

% estimate seasonal indicator
%----------------------------
out.sc=parametric_seasonality(xt);

%  deseasonalize the original series
%------------------------------------
out.sa=y-out.sc;

% irregular component
%---------------------
out.ic=out.sa-out.trend;

out=ts.decomp_format_output(out,y0,doLog);

    function check_inputs()
        
        if ~isa(y0,'ts')||~get(y0,'NumberOfPages')==1
            
             error('first input must be a time series with a single page')
            
        end
        
        if n<2,doLog=[]; end
        
        if isempty(doLog),doLog=false;end
        
        if ~islogical(doLog)
            
            error('second input (doLog) must be a logical')
            
        end
        
        if n<3,dorder=[]; end
        
        if isempty(dorder),dorder=2;end
        
        if ~(isa(dorder,'double') && ...
                isscalar(dorder) && ...
                isfinite(dorder) && ...
                dorder>0 && ...
                ceil(dorder)==floor(dorder))
            
            error('second input (dorder) must be a positive integer')
            
        end
        
    end

    function seas=parametric_seasonality(xt)
        
        freq=utils.time_series.freq2freq(get(y0,'frequency'));
        
        Dum=zeros(T,freq);
        
        for id=1:freq
            
            pos=id:freq:T;
            
            Dum(pos,id)=1;
            
        end
        
        seas=y;
        
        for ii=1:nvar
            
            [~,~,R]=regress(xt(:,ii),Dum);
            
            seas(:,ii)=xt(:,ii)-R;
            
        end
        
        seas=bsxfun(@minus,seas,mean(seas,1));
        
    end

    function trend=parametric_trend()
        
        X=zeros(T,dorder+1);
        
        X(:,1)=1;
        
        for io=1:dorder
            
            X(:,1+io)=(1:T).^io;
            
        end
        
        trend=y;
        
        for ii=1:nvar
            
            [~,~,R]=regress(y(:,ii),X);
            
            trend(:,ii)=y(:,ii)-R;
            
        end
        
    end

end