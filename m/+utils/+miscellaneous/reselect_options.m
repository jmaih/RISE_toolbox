function opt=reselect_options(options,fn)
% reselect_options - discards the options not used by a function
%
% Syntax
% -------
% ::
%
%   opt=reselect_options(options,fn)
%
% Inputs
% -------
%
% - **options** [struct]: all options
%
% - **fn** [char|function_handle]: function whose options we are interested
% in.
%
% Outputs
% --------
%
% - **opt** [struct]: options of function fn
%
% More About
% ------------
%
% - the function of interest should be such that if called without inputs,
% it returns its default options.
%
% Examples
% ---------
%
% See also: 

if ischar(fn)
    fn=str2func(fn);
end

defaults=fn();
fn_fields=fieldnames(defaults);
opt_fields=fieldnames(options);
opt=rmfield(options,opt_fields-fn_fields);

end