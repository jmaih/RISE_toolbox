function varargout=loop_forecast(varargin)
% LOOP_FORECAST - conditional forecast for regime-switching models in a
% loop
%
% More About
% ------------
%
% - This function is the same as RSCF.FORECAST except that it calls
% RSCF.FORECAST in a loop, permitting variations in the unconditional
% regimes.
%
% Examples
% ---------
%
% See also: RSCF.FORECAST

nout=nargout;
if nargin==0
    nf=nargout('utils.forecast.rscond.forecast');
    output=cell(1,nf);
    [output{1:nf}]=utils.forecast.rscond.forecast(varargin{:});
    if nout
        varargout=output(1:nout);
    else
        disp(output{1})
    end
else
    narginchk(5,6)
    
    options=varargin{5};
    varargin{5}.forecast_conditional_sampling_ndraws=1;
    forecast_conditional_sampling_ndraws=options.forecast_conditional_sampling_ndraws;
    for idraw=1:forecast_conditional_sampling_ndraws
        [out{1:nout}]=utils.forecast.rscond.forecast(varargin{:});
        if idraw==1
            varargout=out;
        else
            for iout=1:nout
                varargout{iout}(:,:,idraw)=out{iout};
            end
        end
    end
    for iout=1:nout
        varargout{iout}=squeeze(varargout{iout});
    end
end

end