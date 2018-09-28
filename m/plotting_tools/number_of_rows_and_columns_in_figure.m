function [Remains,r,c]=number_of_rows_and_columns_in_figure(fig,nvar,r0,c0)
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

Remains=nvar-(fig-1)*(r0*c0);

if Remains<r0*c0
    if Remains==3 && r0>=3
        r=3;
        c=1;
    else
        r=min(1,Remains);
        c=min(1,Remains);
        while r*c<Remains
            r=min(r+1,r0);
            if r*c<Remains
                c=min(c+1,c0);
            end
        end
    end
else
    r=r0;
    c=c0;
end