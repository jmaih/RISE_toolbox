function this=exp(this)
% Here it does not make sense to have names any more. But
% all the same, perhaps I should have a function to rename
% the series?
this=rise_time_series(this.start,exp(double(this)));
end
