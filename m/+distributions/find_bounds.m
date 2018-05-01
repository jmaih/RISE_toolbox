function bounds=find_bounds(distr,mm,ss,prob,varargin)
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


[~,~,icdfn,~,m2h]=distributions.(distr)();

[a,b]=m2h(mm,ss,varargin{:});

low=icdfn(.5*(1-prob),a,b,varargin{:});
high=icdfn(.5*(1+prob),a,b,varargin{:});
bounds=[low,high];

end


