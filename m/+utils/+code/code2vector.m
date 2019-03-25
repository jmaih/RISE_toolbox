%  INTERNAL FUNCTION: transforms a set of function handles into a single function
% 
%  ::
% 
%    xout=code2vector(xcell)
% 
%  Args:
% 
%     - **xcell** [fhandle|cell|struct]:
% 
%       - if cell, it is assumed that all entries are anonymous functions
%         sharing the same inputs. One anonymous function is returned in a cell
%       - if fhandle, the same input is returned
%       - if struct
% 
%           - derivatives are identified by the fields
%             'size','functions','map','partitions'
%           - eval items are identified by the fields 'code','argins','argouts'
%           - else it is assumed we have transition matrix, in which case the
%             output is the same as the input
% 
%     - **devect** [true|{false}]: de-vectorize functions.
% 
%  Returns:
%     :
% 
%     - **xout** : vectorized function handle
% 
%  Note:
% 
%     - The routine checks whether the input has ':' or not to writes a
%       consistent unique function
% 
%  Example:
% 
%  See also:
%     -utils.code.code2func
%