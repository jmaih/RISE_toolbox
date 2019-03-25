% SETFIELD Set structure field contents.
%    S = SETFIELD(S,FIELD,V) sets the contents of the specified field to the
%    value V.  For example, SETFIELD(S,'a',V) is equivalent to the syntax
%    S.field = V, and sets the value of field 'a' as V. S must be a 1-by-1
%    structure.  FIELD can be a character vector or string scalar. The
%    changed structure is returned.
% 
%    S = SETFIELD(S,{i,j},'a',{k},V) is equivalent to the syntax
%        S(i,j).field(k) = V;
% 
%    In other words, S = SETFIELD(S,sub1,sub2,...,V) sets the
%    contents of the structure S to V using the subscripts or field
%    references specified in sub1,sub2,etc.  Each set of subscripts in
%    parentheses must be enclosed in a cell array and passed to
%    SETFIELD as a separate input.  Field references are passed as
%    strings or character vectors.  
% 
%    For improved performance, when setting the value of a simple 
%    field, use <a href="matlab:helpview([docroot '/techdoc/matlab_prog/matlab_prog.map'], 'dynamic_field_names')">dynamic field names</a>.
% 
%    See also GETFIELD, ISFIELD, FIELDNAMES, ORDERFIELDS, RMFIELD.
%
%    Reference page in Doc Center
%       doc setfield
%
%    Other functions named setfield
%
%       fints/setfield
%