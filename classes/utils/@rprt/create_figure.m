%--- help for rprt/create_figure ---
%
%  CREATE_FIGURE - Create a figure with multiple subplots.
% 
%  Syntax:
% 
%    obj.create_figure(Subplots)
% 
%    obj.create_figure(Subplots,optName1,optValue1,...)
% 
%  Description:
% 
%    This method creates a figure with multiple subplots based on the given inputs.
%    The subplots can contain different plots with optional titles, legends, plotting
%    functions, data, highlighted dates, and style descriptions.
% 
%  Inputs:
% 
%    - obj : The report object.
% 
%    - varargin : Optional name-value pairs to customize the figure and subplots.
% 
%        - 'DateFormat' : Optional. Date format for the x-axis labels.
% 
%        - 'range' : Optional. Range of values to plot on the x-axis.
% 
%        - 'subplots' : Required. cell array of subplots to add to the figure.
% 
%            Each subplot is defined as a cell within the cell array. The
%            first element of the subplot cell is invariably the data to
%            plot. The remaining elements are optional arguments coming in
%            pairs and defined as follows:
% 
%            - 'caption' : Optional. Title of the subplot.
%            - 'legend' : Optional. legend for the subplot.
%            - 'plotfunc' : Optional. Plotting function to use (default: @plot).
%            - 'highlighted_dates' : Optional. Cell array of 1-by-2 vectors containing
%                the dates to be highlighted on the plot.
%            - 'style' : Optional. Description of the plot style (e.g., linewidth, color, etc.).
%            - 'yline' : Optional. horizontal line on the y axis
%            - 'xline' : Optional. vertical line(s) on the x-axis
% 
%        - other options of rprt.figure : i.e. 'caption', 'width', 'placement',
%          'height', 'angle', 'scale', 'numbering'
% 
%  Example:
%  <<
%    report = rprt('Sample Report', 'John Doe', 'portrait');
%  subplots={
%      {cumsum(randn(100,1)),'caption','the first plot','highlighted_dates',{[20,30]}}
%      {cumsum(randn(100,1)),'caption','the second plot','highlighted_dates',{[20,30]},'style',{'linewidth',2}}
%      };
%    report.createFigure(subplots,'caption', 'My Figure', 'DateFormat', 'dd-mm-yyyy', ...
%        'range', [1, 10]);
%  >>
%