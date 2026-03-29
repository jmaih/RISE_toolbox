%--- help for rsymbdiff.stamper ---
%
%  STAMPER Creates stamping, querying, and assigning functions for managing derivatives.
% 
%    [stamp, query, assign, status] = STAMPER(TMP) returns three functions
%    and a containers.Map object that together facilitate the management of symbolic
%    derivatives.
% 
%    Input arguments:
%    - TMP: prefix for temporary terms
% 
%    Output arguments:
%    - stamp: A function for generating unique stamping_engines in the format "T_A," where A
%             is an incrementing integer.
%    - query: A function for querying the existence of a stamping_engine (e.g., "T_A_B") and
%             retrieving its associated value. It keeps track of how many times
%             a stamping_engine has been queried.
%    - assign: A function for assigning a value to a stamping_engine (e.g., "T_A_B") in a
%              containers.Map object. It also resets the query_engine count for that stamping_engine.
%    - status: A function that summarizes the status of the computations.
% 
%    Usage:
%    [stamp, query, assign, status] = stamper();
% 
%    Example:
%    [stamp, query, assign, status] = stamper();
%    stamping_1 = stamp();       % Generates a new stamping_engine "T_1"
%    assign('T_1_B', 'T_1');    % Assigns "T_1" to "T_1_B" in the map
%    value = query('T_1_B');     % Queries the value associated with "T_1_B"
%    computation_status = status();  % Retrieves computation status
% 
%    See also: containers.Map
%