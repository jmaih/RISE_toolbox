%  APPEND Combine string elements
%     STR = APPEND(STR1,STR2) combines elements of STR1 and STR2 into STR.
%     Each input argument can be a string, character vector, or a cell array
%     of character vectors.
%     * If any input is a string array, then STR is a string array.
%     * If no input is a string array, and at least one is a cell array of
%       character vectors, then STR is a cell array of character vectors.
%     * If all inputs are character vectors, then STR is a character vector.
% 
%     STR1 and STR2 must have compatible sizes. Two inputs have compatible
%     sizes if, for every dimension, the dimension sizes of the inputs are
%     either the same or one of them is 1. In the simplest cases, they can
%     be the same size or one can be a scalar.
%  
%     STR = APPEND(STR1,STR2,STR3,...) combines elements of all inputs using
%     the rules specified in the previous syntax.
% 
%     Example:
%         append("data", ".tar.gz")
%  
%         returns  
%  
%             "data.tar.gz"
%   
% 
%         append(["paper1","paper2"], '.docx')
%   
%         returns  
%  
%             "paper1.docx"    "paper2.docx"
%   
% 
%         % Combine a char, column cell, and row string
%         append('data', {'.dat';'.tar'}, ["",".gz"])
%  
%         returns
%  
%             "data.dat"    "data.dat.gz"
%             "data.tar"    "data.tar.gz"
%  
%     See also STRING/PLUS, STRCAT, HORZCAT, JOIN
%
%    Documentation for append
%       doc append
%
%    Other uses of append
%
%       matlab.net.http.io.ContentConsumer/append
%       matlab.unittest.internal.diagnostics.CompositeConditionsSupplier/append
%       matlab.unittest.plugins.plugindata.ResultDetails/append
%       rprt/append
%       timeseries/append
%       TreeBagger/append
%