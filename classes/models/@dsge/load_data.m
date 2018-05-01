function varargout=load_data(obj)
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

just_starting=isempty(obj);

[varargout{1:nargout}]=load_data@generic(obj);

if just_starting
    
    return
    
end

obj=varargout{1};
data_provided=obj.options.data.NumberOfVariables>0;
if data_provided %simulation_available || 
    estim_start_date=obj.options.estim_start_date;
    obj.dates_filtering=date2serial(estim_start_date);
    obj.dates_smoothing=obj.dates_filtering;
%    obj.dates_filtering=serial2date(date2serial(estim_start_date)+(0:obj.data.nobs));
%    obj.dates_smoothing=obj.dates_filtering(1:end-1);
     
    % information on the conditional variables
    %------------------------------------------
    kmax=max(obj.exogenous.shock_horizon(:));
    if kmax
        if obj.data.npages<kmax+2
            warning('the anticipation horizon of agents will be reduced since it exceeds the number of advance information')
        end
        obj.exogenous.shock_horizon=min(obj.exogenous.shock_horizon,kmax+2);
    end
end

varargout{1}=obj;

end