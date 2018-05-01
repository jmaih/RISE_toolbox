function [mvcb,gam_,my]=chebyshev_box(this,gam)
% chebyshev_box - constructs chebyshev boxes for multivariate-multiperiods
% densities
%
% ::
%
%
%   mvcb=chebyshev_box(this,gam)
%
% Args:
%
%    - **this** [ts|rts] : time series with many pages (number of simulations)
%    and potentially many columns (number of variables)
%
%    - **gam** [scalar|vector] : percentile(s)
%
% Returns:
%    :
%
%    - **mvcb** [struct] : structure with time series in its fields. Each
%    field represents a particular variable
%
%    - **gam_** [scalar|vector] : sorted percentile(s)
%
%    - **my** [ts] : mean across simulations
%
% Note:
%
% Example:
%
%    See also:

Y=double(this);
Y=permute(Y,[3,1,2]);
gam_=sort(gam);
[mvcb,my]=utils.forecast.kolsrud.multivariate_chebyshev_box(Y,gam_);
% fourth dimension is the number of gam(s)
mvcb_=permute(mvcb,[2,3,1,4]);
mvcb=struct();
varnames=get(this,'varnames');
start=get(this,'start');
G=get(this,'NumberOfVariables');
my=ts(start,permute(my,[2,3,1]),varnames);
for g=1:G
    newdata=squeeze(mvcb_(:,g,:,:));
    if g==1
        prototype=ts(start,newdata,{'low','high'});
    end
    mvcb.(varnames{g})=reset_data(prototype,newdata);
end

end
