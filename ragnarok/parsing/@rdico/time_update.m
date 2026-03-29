%--- help for rdico.time_update ---
%
%  construct a function handle that updates time
% 
%  Inputs:
% 
%    var_list: the list of eligible variables
% 
%    default_time: {t|''}
% 
%  Output:
% 
%    updater: function handle that can be called as updater(expr,inc_flag)
%    where :
% 
%    - expr is the expression to update.
%    - inc_flag : true for increments, false for decrements
%