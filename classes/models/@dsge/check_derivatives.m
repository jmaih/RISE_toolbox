%--- help for dsge/check_derivatives ---
%
%  Compares the derivatives and the solutions from various differentiation techniques
% 
%  ::
% 
%    check_derivatives(obj)
%    retcode=check_derivatives(obj)
% 
%  Args:
%     obj (rise | dsge): model object or vectors of model objects
% 
%  Returns:
%     :
% 
%     - **retcode** [numeric]: 0 if no problem is encountered during the
%       comparisons. Else the meaning of recode can be found by running
%       decipher(retcode)
% 
%  Note:
% 
%     - The derivatives computed are 'automatic', 'symbolic' or 'numerical'
%     - The comparisons are done relative to automatic derivatives, which are
%       assumed to be the most accurate.
% 
%