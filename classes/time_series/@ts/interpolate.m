function this=interpolate(this,method,varargin)
% Interpolate based on the time series to fill in the missing dates
%
% ::
%
%    db = interpolate(db);
%    db = interpolate(db, method, varargin);
%
% Args:
%
%    db (ts object): time series object to interpolate from
%
%    method (string): interpolate method (default: spline). Same as the option used in MATLAB, i.e.,
%
%        - 'nearest': nearest neighbor interpolation
%        - 'linear': linear interpolation
%        - 'spline': piecewise cubic spline interpolation (SPLINE)
%        - 'pchip': shape-preserving piecewise cubic interpolation
%        - 'cubic': same as 'pchip'
%        - 'v5cubic': the cubic interpolation from MATLAB 5, which does not
%          extrapolate and uses 'spline' if X is not equally spaced.
%
%    varargin{1}: optional condition for extrapolation (default: no extrapolation)
%
%       - 'extrap': extrapolate out for dates outside the dates with observations
%       - extrapval (double): use extrapval for dates outside the dates with observations
%
    
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
