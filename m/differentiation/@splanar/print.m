%    PRINT Print or save a figure or model.
%      A subset of the available options is presented below. For more details
%      see <a href="matlab:helpview([docroot '/matlab/ref/print.html'])" />the documentation</a>.
% 
%      PRINT, by itself, prints the current figure to your default printer.
%      Use the -s option to print the current model instead of the current figure.
%        print         % print the current figure to the default printer
%        print -s      % print the current model to the default printer
% 
%      PRINT(filename, formattype) saves the current figure to a file in the
%      specified format. Vector graphics, such as PDF ('-dpdf'), and encapsulated
%      PostScript ('-depsc'), as well as images such as JPEG ('-djpeg') and PNG ('-dpng')
%      can be created. Use '-d' to specify the formattype option
%        print(fig, '-dpdf', 'myfigure.pdf'); % save to the 'myfigure.pdf' file
%      The full list of formats is <a href="matlab:helpview([docroot '/matlab/ref/print.html#inputarg_formattype'])" />documented here</a>.
% 
%      PRINT(printer, ...) prints the figure or model to the specified printer.
%      Use '-P' to specify the printer option.
%        print(fig, '-Pmyprinter'); % print to the printer named 'myprinter'
% 
%      PRINT(resize,...) resizes the figure to fit the page when printing.
%      The resize options are valid only for figures, and only for page
%      formats (PDF, and PS) and printers. Specify resize as either
%        '-bestfit'  to preserve the figure's aspect ratio or
%        '-fillpage' to ignore the aspect ratio.
% 
%    <a href="matlab:helpview([docroot '/matlab/ref/print.html'])" />The documentation</a> contains additonal details and examples, including how to
%    specify the figure or model to print, adjust the output size and
%    resolution, save to the clipboard, and specify the renderer to use.
% 
%    See also SAVEAS, PRINTPREVIEW, SAVEFIG.
%
%    Reference page in Doc Center
%       doc print
%
%    Other functions named print
%
%       arima/print    splanar/print
%