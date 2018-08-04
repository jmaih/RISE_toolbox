function [trend,detrended]=ma_filter(varargin)
% MA_FILTER - moving average filter
%
% ::
%
%
%   [trend,detrended]=ma_filter(y,q)
%   [trend,detrended]=ma_filter(y,q,extend)
%
% Args:
%
%    - **y** [ts] : scalar time series
%      the time series. The following are admitted:
%      - annual data : e.g. 1990 or '1990'
%      - bi-annual data : e.g. '1990H1'
%      - Quarterly data : e.g. '1990Q3'
%      - monthly data : e.g. '1990M12'
%    - **q** [integer|{0.5*frequency}] : number of periods before or after the
%      current one to be considered in the moving average calculation. The
%      total window length is 2q+1
%    - **extend** [true|{'false'}] : if true, replicated observations are
%      added both at the beginning and at the end of the original dataset in
%      order to avoid losing some observations during the filtering process.
%
% Returns:
%    :
%
%    - **trend** [ts] : (non-parametric) trend
%
%    - **detrended** [ts] : y-trend
%
% Note:
%
% Example:
%
%    See also:
%    ---------
%

n=length(varargin);

y=varargin{1};

q=[];

extend=[];

check_inputs()

varargin{1}=double(y);

trend=utils.time_series.mymafilter(double(y),q,extend);

if extend
    
    trend=ts(y.start_date_number,trend);
    
else
    
    trend=ts(y.start_date_number+q,trend);
    
end

if nargout>1
    
    detrended=y-trend;
    
end

    function check_inputs()
        
        if ~isa(y,'ts')
            
            error('first input must be a scalar ts object')
        
        end
        
        if n>1,q=varargin{2}; end
        
        if isempty(q)
            
            qtmp=utils.time_series.freq2freq(get(y,'frequency'));
            
            q=qtmp*0.5;
            
        end
        
        if n>2,extend=varargin{3}; end
        
        if isempty(extend),extend=false; end
        
        if ~isa(extend,'logical')
            
            error('third input (extend) must be true or false')
            
        end
                
        if n>3
            
            error('too many input arguments')
            
        end
        
    end

end