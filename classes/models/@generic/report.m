%--- help for generic/report ---
%
%  Assigns the elements of interest to a rise_report.report object
% 
%    - display the default inputs::
% 
%         report(rise.empty(0))
% 
%    - assign the reported elements in rep_items to destination_root::
% 
%         report(obj,destination_root,rep_items)
% 
%    - assign varargin to obj before doing the rest::
% 
%         report(obj,destination_root,rep_items,varargin)
% 
%  Args:
% 
%     obj (rise | dsge):
%     destination_root (rise_report.report): handle for the actual report
%     rep_items (char | cellstr): list of desired items to report. This list
%        can only include
% 
%        - endogenous
%        - exogenous
%        - observables
%        - parameters
%        - solution
%        - estimation
%        - estimation_statistics
%        - equations
%        - code
% 
%  Returns:
%     :
% 
%     - none (Update the report in place)
% 
%