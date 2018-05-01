function [obs,flag_]=date2obs(start_date,new_date)
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

[start_serial,start_frequency]=date2serial(start_date);

[new_serial,new_frequency]=date2serial(new_date);

if ~strcmpi(start_frequency,new_frequency)
    
    error('dates are not in the same frequency')
    
end

retcode=1-(new_serial>=start_serial);

if nargout>1
    
    flag_=retcode;
    
end

if retcode
    
    if nargout<2
        
        error('new date occurs before the start date')
        
    end
    
    obs=[];
    
else
    
    obs=nan(size(new_serial));
    
    for iserial=1:numel(obs)
        
        span=start_serial:new_serial(iserial);
        
        obs(iserial)=numel(span);
        
    end
    
end

end

%     start=The_data.TimeInfo(1).date_2_observation(obj.options.estim_start_date);
%     finish=The_data.TimeInfo(1).date_2_observation(obj.options.estim_end_date);
