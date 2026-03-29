%--- help for dsge/plan2database ---
%
%  plan2database creates a database suitable for conditional forecasting
% 
%  ::
% 
%    db=plan2database(obj,plan,start_date,end_date)
% 
%  Args:
% 
%     obj (rise | dsge): model file
% 
%     *plan* : n x 3 cell array where n is the number of conditions
%      - the first column includes the names of the variables in strings or in
%        cell array
%      - the second column includes the dates for conditioning
%      - the third column includes the values scalar or vectors
% 
%     *start_date*: start date of the database as recognizable by the ts
%        class. this is typically the end of history
% 
%     *end_date*: end date of the database as recognizable by the ts
%        class.
% 
% 
%  Returns:
%     :
% 
%     - **db** [ts]: time series object with the various variables in and
%       conditioning information
% 
%  N.B:
%     :
% 
%       There cannot be multiple variables and multiple dates
%       simultaneously in one row
% 
%  See also ts.fold, dsge.growth_database, dsge.initial_conditions
%