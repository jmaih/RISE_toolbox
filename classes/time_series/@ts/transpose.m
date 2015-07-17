function d=transpose(db)
% H1 line
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
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

d=main_frame(db,false);

d=permute(d,[2,1,3]);
end