function [r,c]=one_figure_rows_columns(nitems)
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

r=ceil(sqrt(nitems));
c=floor(sqrt(nitems));
odd=true;
while r*c<nitems
    if odd
        c=c+1;
    else
        r=r+1;
    end
    odd=~odd;
end

end
