%--- help for dsge/is_stationary_system ---
%
%  Checks whether a the model is stationary. i.e., does not containt trending variables
% 
%  ::
% 
%    flag=is_stationary_system(obj)
% 
%  Args:
% 
%     obj (rise|dsge): rise/dsge model object
% 
%  Returns:
%     :
% 
%     - **flag** (bool): true if the model is stationary
% 
%  Note:
% 
%     - There is a difference between stability and stationarity
% 
%        - stability refers to the system as a whole and conditions for
%          stability are often assessed through eigenvalues inside the unit circle
%        - stationarity refers to a scalar stochastic process. And such a
%          process will be said stationary if its first and second moment do not
%          vary with time.
% 
%  See also:
% 
%     - rise_generic/is_stable_system
% 
%