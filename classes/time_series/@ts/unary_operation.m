function db=unary_operation(db,op_string)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

if ischar(op_string)
    op_string=str2func(['@(x)',op_string,'(x)']);
elseif ~isa(op_string,'function_handle')
    error('function must be a string or a function handle')
end
db=ts(db.date_numbers,op_string(db.data));
end
