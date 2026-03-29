% RECESSION_DATES Detects recession dates based on the given frequency.
%    r = recession_dates(F)
% 
%    Syntax:
%        r = recession_dates(F)
% 
%    Description:
%        This function detects recession dates based on the given frequency.
%        The function retrieves recession data from NBER (National Bureau of
%        Economic Research) and returns a cell array containing the start
%        and end dates of each recession period.
% 
%    Inputs:
%        - F: Frequency of recession detection, specified as a character.
%            Valid options are:
%                'D' - Daily
%                'W' - Weekly
%                'M' - Monthly
%                'Q' - Quarterly
%                'H' - Half-yearly
%                'Y' - Yearly
%            If F is not provided, the default frequency is 'Q' (Quarterly).
% 
%    Outputs:
%        - r: Recession dates, returned as a cell array of size Nx2, where N
%            is the number of recession periods detected. Each row of the
%            cell array represents a recession period and contains the start
%            date in the first column and the end date in the second column.
% 
%    Examples:
%        % Detect quarterly recession dates
%        r1 = recession_dates('Q');
% 
%        % Detect monthly recession dates
%        r2 = recession_dates('M');
% 
%        % Detect yearly recession dates
%        r3 = recession_dates('Y');
% 
%    See also: fetch_fred, convert
%