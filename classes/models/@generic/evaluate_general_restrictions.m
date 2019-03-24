%--- help for generic/evaluate_general_restrictions ---
%
%  Evaluate general restrictions
% 
%  ::
% 
%    g = evaluate_general_restrictions(obj)
% 
%  Args:
% 
%     obj (rise | dsge | rfvar | svar): scalar of vector or RISE model objects
% 
%  Returns:
%     :
% 
%     - **g** [cell array] : value of restrictions for each object.
% 
%  Note:
% 
%     - The restrictions will be processed as g(x)<=0. But all the user has to
%       do is to put zero where the restrictions are not violated!!!
% 
%