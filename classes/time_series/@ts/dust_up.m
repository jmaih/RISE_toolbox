function this=dust_up(this,crit)
% dust_up - sets insignificant digits to zero
%
% ::
%
%
%   this=dust_up(this)
%
%   this=dust_up(this,crit)
%
% Args:
%
%    - **this** [rts|ts]: time-series object
%
%    - **crit** [numeric|{sqrt(eps)}]: cutoff point
%
% Returns:
%    :
%
%    - **this** [rts|ts]: time-series object
%
% Note:
%
% Example:
%
%    See also:

if nargin<2
    crit=sqrt(eps);
end

this.data(abs(this.data)<crit)=0;

end
