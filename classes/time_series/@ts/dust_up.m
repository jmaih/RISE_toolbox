function this=dust_up(this,crit)
if nargin<2
    crit=sqrt(eps);
end

% dust_up - sets insignificant digits to zero
%
% Syntax
% -------
% ::
%
%   this=dust_up(this)
%
%   this=dust_up(this,crit)
%
% Inputs
% -------
%
% - **this** [rts|ts]: time-series object
%
% - **crit** [numeric|{sqrt(eps)}]: cutoff point
%
% Outputs
% --------
%
% - **this** [rts|ts]: time-series object
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 
this.data(abs(this.data)<crit)=0;

end
