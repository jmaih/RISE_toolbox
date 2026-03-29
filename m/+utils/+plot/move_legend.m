%  MOVE_LEGEND repositions a legend within a figure.
% 
%  Syntax:
%    move_legend(legdHandle, nrows, ncols, w)
% 
%  Inputs:
%    - legdHandle: Handle to the legend that you want to reposition.
%    - nrows: Number of rows in the subplot grid.
%    - ncols: Number of columns in the subplot grid.
%    - w: Subplot position where the legend should be centered.
% 
%  Description:
%    This function allows you to move the position of a legend within a figure
%    containing a subplot grid defined by the number of rows and columns. You
%    specify the subplot position (w) where the legend should be centered.
% 
%  Example:
%    % Create a sample figure with subplots
%    figure;
%    subplot(2, 2, 1);
%    plot(1:10, rand(1, 10), 'r');
%    hold on;
%    plot(1:10, rand(1, 10), 'b');
%    legend('Red', 'Blue');
% 
%    % Move the legend to the center of subplot(2, 2, 3)
%    move_legend(legend, 2, 2, 3);
% 
%  See also:
%    legend, subplot
%