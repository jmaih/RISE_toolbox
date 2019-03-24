%--- help for load_data ---
%
%  INTERNAL FUNCTION: loads data for estimation and for forecasting
% 
%  ::
% 
%     [obj,issue,retcode] = load_data(obj)
%     [obj,issue,retcode] = load_data(obj,varargin)
% 
%  Args:
% 
%     obj (rise | dsge | rfvar | svar): scalar or vector or RISE model objects
%     data (ts | struct): data to load
%     data_demean (true | {false}): flag for demeaning the data
% 
%  Returns:
%     :
% 
%     - **obj** [rise|dsge|rfvar|svar]: scalar or vector or RISE model objects
%     - **issue** [{''}|char]: message describing any issue related to the
%       loading of the data
%     - **retcode** [{0}|positive integer]: if retcode is different from 0,
%       then there is a problem
%
%    Other functions named load_data
%
%       dsge/load_data    generic/load_data
%