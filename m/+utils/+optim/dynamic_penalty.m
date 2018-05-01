function pp=dynamic_penalty(vector,k,alpha_)
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

if nargin<3
    alpha_=100;
end
xi=max(0,vector);
pp=d_(k)*H_(xi);
    function h=H_(xi)
        h=sum(theta_(xi).*xi.^gamma_(xi));
        function tt=theta_(xi)
            tt=xi;
            tt(xi<1e-5)=alpha_;
            tt(xi>=1e-5 & xi<1e-3)=10*alpha_;
            tt(xi>=1e-3 & xi<1)=100*alpha_;
            tt(xi>=1)=1000*alpha_;
        end
        function gg=gamma_(xi)
            gg=2*ones(size(xi));
            gg(xi<1)=1;
        end
    end
    function d=d_(k)
        d=k^0.1;
    end
end
