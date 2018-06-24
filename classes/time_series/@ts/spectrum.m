function [sw,jj,T]=spectrum(this)
% Computes the spectral density of the data
%
% ::
%
%    [sw,jj,T]=SPECTRUM(this)
%
% Args:
%
%    - **this** [ts]: time series object
%
% Returns:
%
%    - **sw** [matrix|vector]: spectrum of potentially multiple time series
%      in columns
%    - **jj** [vector]: range of the spectrum (x-axis)
%    - **T** [integer]: number of observations
%

y=this.data;

if size(y,3)>1

    error([mfilename,':: this operation is only defined for databases with one page'])

end

bad=any(isnan(y),2);

y=y(~bad,:);

if isempty(y)

    error([mfilename,':: no valid data to construct the periodogram'])

end

[T,n]=size(y);

% autocovariances
y=bsxfun(@minus,y,mean(y,1));

gam=nan(T,n);

for ii=1:T

    gam(ii,:)=sum(y(ii:T,:).*y(1:T-ii+1,:),1);

end

gam=gam/T;

% spectrum
jj=(1:T-1)';

w=2*pi*jj/T;

sw=nan(T-1,n);

for ii=1:numel(w)

    cwj=cos(w(ii)*jj);

    sw(ii,:)=1/(2*pi)*(gam(1,:)+2*sum(bsxfun(@times,gam(2:end,:),cwj),1));

end

end
