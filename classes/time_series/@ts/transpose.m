function d=transpose(db)
d% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

=main_frame(db,false);

d=permute(d,[2,1,3]);
end