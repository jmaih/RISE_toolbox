function vcov=farthest(vcov0,debug)
% INTERNAL FUNCTION
%

if nargin<2
    debug=false;
end

vcov=utils.cov.nearest(vcov0,debug,true);

end