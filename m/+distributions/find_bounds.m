function bounds=find_bounds(distr,mm,ss,prob,varargin)
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


[~,~,icdfn,~,m2h]=distributions.(distr)();

[a,b]=m2h(mm,ss,varargin{:});

low=icdfn(.5*(1-prob),a,b,varargin{:});
high=icdfn(.5*(1+prob),a,b,varargin{:});
bounds=[low,high];

end


