function [Remains,r,c]=number_of_rows_and_columns_in_figure(fig,nvar,r0,c0)
Remains=nvar-(fig-1)*(r0*c0);

if Remains<r0*c0
    r=min(1,Remains);
    c=min(1,Remains);
    while r*c<Remains
        r=min(r+1,r0);
        if r*c<Remains
            c=min(c+1,c0);
        end
    end
else
    r=r0;
    c=c0;
end