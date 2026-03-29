%  OPTIMZE_LOGISTIC -- finds the optimal parameters of a logistic function
% 
%  ::
% 
% 
%    ab=optimize_logistic(x1,p1,x2,p2)
% 
%    ab=optimize_logistic(x1,p1,x2,p2,op)
% 
%  Args:
% 
%     - **x1,x2** [scalars] : points at which to evaluate the logistic
% 
%     - **p1,p2** [scalars] : probability values for x1 and x2
% 
%     - **op** [{'+'}|'-'] : decides whether b is to be multiplied by 1 or by
%       -1 in the evaluation the logistic function.
% 
%  Returns:
%     :
% 
%     - **ab** [vector] : optimized parameters of the logistic function
% 
%  Note:
% 
%       The function has the form f(x,g,c)=1/(1+exp(-g*(x-c))) so that
%       f(x1,g,c)=p1 and f(x2,g,c)=p2
% 
%  Example:
% 
%     See also:
%