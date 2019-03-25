%  multiple plots in graphs
% 
%  ::
% 
%    hfig=multiple(plotfunc,vnames,fig_title,r0,c0)
%    hfig=multiple(plotfunc,vnames,fig_title,r0,c0,varargin)
% 
%  Args:
% 
%     - **plotfunc** [function handle]: which takes as input a valid variable
%       name and returns (1) the name/description of the variable/parameter
%       plotted, (2) the legend
%     - **vnames** [cellstr]: names of variables to be plotted
%     - **fig_title** [char]: main title of the figures to plot
%     - **r0** [integer]:: the desired maximum number of rows in each figure
%     - **c0** [integer]:: the desired maximum number of columns in each figure
%     - **varargin** [|{}]:: pairwise elements of entering the title and the
%       legend output arguments
% 
%  Returns:
%     :
% 
%     - **hfig** [handles]: handles to the different figures created
% 
%