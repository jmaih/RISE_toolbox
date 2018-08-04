function this=interpolate(this,method,varargin)
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


%     dbq = interpolate(db) interpolates to find dbq, the values of the
%     underlying function V=F(X) at the query points found to be the
%     missing values in the time series, which can be univariate or
%     multivariate and can have many pages
%  
%     dbq = interpolate(db,METHOD) specifies alternate methods.
%     The default is spline interpolation. Available methods are:
%  
%       'nearest'  - nearest neighbor interpolation
%       'linear'   - linear interpolation
%       'spline'   - piecewise cubic spline interpolation (SPLINE)
%       'pchip'    - shape-preserving piecewise cubic interpolation
%       'cubic'    - same as 'pchip'
%       'v5cubic'  - the cubic interpolation from MATLAB 5, which does not
%                    extrapolate and uses 'spline' if X is not equally
%                    spaced.
%  
%     dbq = interpolate(db,METHOD,'extrap') uses the interpolation algorithm
%     specified by METHOD to perform extrapolation for the missing values
%     for which the dates are outside the dates of the non-missing values.
%  
%     dbq = interpolate(db,METHOD,EXTRAPVAL) replaces the values outside of
%     the interval spanned by X with EXTRAPVAL.  
    
if nargin<2 || isempty(method)
    method='spline'; 
end

dn=this.date_numbers;
v=this.data;
[~,nvar,npages]=size(v);
for ivar=1:nvar
    for ipage=1:npages
        query_points=isnan(v(:,ivar,ipage));
        if any(query_points)
            v(query_points,ivar,ipage)=interpolate(dn(~query_points),...
                v(~query_points,ivar,ipage),dn(query_points),method,varargin{:});
        end
    end
end
this.data=v;
end
