function [ishigami,bounds]=ishigami()

ishigami=@my_ishigami;
%%%========================================
%     1		      0.3307
%     2		      0.3309
%     3		      0.0007
%     1     2      0.0217
%     1     3      0.2782
%     2     3      0.0420
%    1.0042
%%%========================================

bounds=repmat([-pi,pi],3,1);

end

function f=my_ishigami(x,debug)
if nargin<2
    debug=false;
end
f=sin(x(1,:))+7*sin(x(2,:)).^2+0.1*x(3,:).^4.*sin(x(1,:));
if debug
    f=f(ones(10,1),:);
end
end