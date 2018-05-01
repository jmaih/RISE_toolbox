function [st,Q0,PAI,retcode]=choose_state(st,Qfunc,PAI,y)
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

% st: state
% Qfunc: function that returns the transition matrix evaluated at y ...
% PAI: current updated probabilities
% y: data for current period
% Q0: is nan if st is given
retcode=0;
Q0=nan(size(PAI,1));
if isempty(st)||isnan(st)
    [Q0,retcode]=Qfunc(y);
    % update probabilities
    %---------------------
    if retcode
        return
    end
    PAI=Q0'*PAI;
    csp=[0;cumsum(PAI)];
    st=find(csp>rand,1,'first')-1;
end
end
