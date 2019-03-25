%  INTERNAL FUNCTION: load the start values for forecasting
% 
%  ::
% 
%    v=load_start_values(names,db,date)
%    v=load_start_values(names,db,date,v)
% 
%  Args:
% 
%     - **names** [char|cellstr]: names of variables to load from the database
%     - **db** [ts|struct]: ts time series data
%     - **date** [char|numeric|serial date]: date for which the data are
%       requested
%     - **v** (optional)[vector,3-dimensional array|[]]: start values. If
%       provided, it should have the same number of rows as the number of names
% 
%  Returns:
%     :
% 
%     - **v** [vector,3-dimensional array|[]]: start values
% 
%  Note:
% 
%     - leads are ignored
%     - if a variable is not found in the database, it is initialized at 0 in
%       the first page and nan in the subsequent pages
% 
%  Example:
% 
%     db=ts('1990q1',rand(100,5,3),{'v1','v2','v3','v4','v5'})
%     v=utils.forecast.load_start_values({'v4','v5','v5_AUX_L_8{+2}','v5_AUX_F_3'},db,'2011Q4')
% 
%  See also:
%     data_request
% 
%