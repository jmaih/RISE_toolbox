% SAVEAS Save Figure or Simulink block diagram in desired output format
%    SAVEAS(H,'FILENAME')
%    Will save the Figure or Simulink block diagram with handle H to file 
%    called FILENAME. 
%    The format of the file is determined from the extension of FILENAME.
% 
%    SAVEAS(H,'FILENAME','FORMAT')
%    Will save the Figure or Simulink block diagram  with handle H to file 
%    called FILENAME in the format specified by FORMAT. FORMAT can be the 
%    same values as extensions of FILENAME. 
%    The FILENAME extension does not have to be the same as FORMAT.  
%    The specified FORMAT overrides FILENAME extension.
% 
%    Valid options for FORMAT are:
% 
%    'fig'  - save figure to a single binary FIG-file.  Reload using OPEN. 
%    'm'    - save figure to binary FIG-file, and produce callable
%             MATLAB code file for reload.
%    'mfig' - same as M.
%    'mmat' - save figure to callable MATLAB code file as series of creation
%             commands with param-value pair arguments.  Large data is saved
%             to MAT-file.  
%             Note: MMAT Does not support some newer graphics features. Use
%                   this format only when code inspection is the primary goal.
%                   FIG-files support all features, and load more quickly. 
% 
%    Additional FORMAT options include devices allowed by PRINT.
% 
%    NOTE: not all format options are allowed for Simulink block diagrams.
%    See the online help for more information.
% 
%    Examples:
% 
%    Write current figure to MATLAB fig file
% 
%        saveas(gcf, 'output', 'fig')
% 
%    Write current figure to windows bitmap file
% 
%        saveas(gcf, 'output', 'bmp')
% 
%    Write block diagram named 'demo' to an Encapsulated Postscript file
% 
%        saveas(get_param('demo', 'Handle'), 'output', 'eps')
% 
%    In the following list, SAVE_SYSTEM is available for Simulink users only. 
%    See also LOAD, SAVE, OPEN, PRINT, SAVE_SYSTEM, EXPORTGRAPHICS.
%
%    Documentation for saveas
%       doc saveas
%
%    Other uses of saveas
%
%       rdico/saveas
%