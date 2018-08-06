function this=sinh(this)
% Overloaded sinh function for ts object
%

% Here it does not make sense to have names any more. But
% all the same, perhaps I should have a function to rename
% the series?
this=ts.unary_operation(this,mfilename);
end