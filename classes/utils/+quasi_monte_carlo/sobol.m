function theta_ = sobol(lb,ub,K,scramble_flag)
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


%   Examples:
%
%      Create a 5-dimensional point set and get the first 1024 points:
%         P = sobolset(5);
%         X = net(P,1024)
%
%      Create a point set and get the 1st, 3rd, 5th... points:
%         P = sobolset(5);
%         X = P([1 3 5 7 9 11 13], :)
%
%      Create a scrambled point set:
%         P = sobolset(5);
%         P = scramble(P,'MatousekAffineOwen');
%         X = net(P,1024)
%
%   See also NET, QRANDSTREAM, SCRAMBLE, SUBSREF.

if nargin<4
    scramble_flag=true;
end

n=size(lb,1);

if exist([mfilename,'set.m'],'file')
    P = sobolset(n);
    if scramble_flag
        P = scramble(P,'MatousekAffineOwen');
    end
    Seq = transpose(net(P,K));
    theta_=bsxfun(@plus,lb,bsxfun(@times,ub-lb,Seq));
else
    theta_ = quasi_monte_carlo.latin_hypercube(n,K,lb,ub);
end
end
