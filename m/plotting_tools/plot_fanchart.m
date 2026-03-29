%  PLOT_FANCHART plots fancharts
% 
%  ::
% 
%     plot_fanchart(data)
%     plot_fanchart(data, MainColor)
%     plot_fanchart(data, MainColor, nticks)
%     hh = plot_fanchart(...)
% 
%  Args:
% 
%     data (struct): Output of ts.fanchart
% 
%     MainColor (char|vector|{'nb'}): Main color for the fanchart.
% 
%        - if char, MainColor must be an element of {'c','b','g','r','m','k','w'}.
%          If the default ('nb') is used, then the number of elements to plot should be exactly 4.
%        - if vector, MainColor must be of the format rgb.
% 
%     nticks (numeric|{8}): Number of tick points in the x-axis.
% 
%  Returns:
% 
%     hh (numeric): Handle to the plot
% 
%  Example:
% 
%     this = ts('1990Q2', rand(100, 1000));
%     out = fanchart(this, [30, 50, 70, 90]);
%     out = fanchart(this, [30, 50, 70, 90]/100);
%     plot_fanchart(out, 'r', 10)
% 
%  See also: FANCHART
%