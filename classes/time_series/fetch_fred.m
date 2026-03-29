%  fetch_fred collects data from the Federal Reserve Economic
%  Data (FRED) database. It allows you to retrieve multiple economic time
%  series data at once. The function accepts a cell array of FRED series IDs
%  and returns a structure containing the retrieved data for each series. 
% 
%  Syntax ::
% 
%    db = fetch_fred(seriesIDs)
% 
%  Inputs:
% 
%  - seriesIDs (cell array of strings): A cell array containing FRED series
%    IDs for the data to be retrieved. Series IDs are unique identifiers for
%    specific economic time series data available on the FRED database. 
% 
%  Outputs:
% 
%  - db (structure array): A structure array containing the retrieved data
%    for each series. The structure has the following fields: 
% 
%    * header (string): The header information from the FRED data file.
%    * Title (string): The title of the time series data.
%    * Frequency (string): The frequency of the data (e.g., Daily, Weekly,
%      Monthly, Quarterly, Annual). 
%    * series (timeseries object): A time series object representing the
%      data, with appropriate time values and corresponding data values. 
% 
%  Examples
%  Example 1: Fetching Single Series
% 
%  seriesID = 'GDPC1'; % Gross Domestic Product (GDP) series ID
%  data = fetch_fred(seriesID);
% 
%  Example 2: Fetching Multiple Series
% 
%  seriesIDs = {'GDPC1', 'UNRATE', 'CPIAUCSL'};
%  data = fetch_fred(seriesIDs);
%