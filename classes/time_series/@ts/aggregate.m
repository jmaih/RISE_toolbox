%--- help for ts/aggregate ---
%
%  Convert time series to a lower frequency
% 
%  ::
% 
%     this=aggregate(this,newfreq);
%     this=aggregate(this,newfreq,method);
% 
%  Args:
%     this (ts object): time series object
%     newfreq (char): frequency to convert the data to. It must be a lower than the frequency of the data in **this**. Available frequencies are
% 
%        - 'Q': quarterly
%        - 'M': monthly
%        - 'H': semiannual
%        - 'W': weekly
%        - 'D': daily
% 
%     method ('interpolation' , 'distribution'): method of conversion (default: 'distribution')
% 
%        - 'interpolation': Use the value corresponding to the first date of the coarse frequency
%        - 'distribution': Use the average of all observations in the sampling frequency
% 
%  Returns:
%     :
% 
%     - **this** (ts object): time series object
% 
%