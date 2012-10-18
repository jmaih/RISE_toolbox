function val=test_system_fun(x,c,optimal_policy)

error(nargchk(1,3,nargin))

nvar=numel(x);
if nargin<3
    optimal_policy=true;
    if nargin<2
        c=1:nvar;
    end
end


y1=c(1)*x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8)+x(9)+x(10);
y2=x(1)+c(2)*x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8)+x(9)+x(10);
y3=x(1)+x(2)+c(3)*x(3)+x(4)+x(5)+x(6)+x(7)+x(8)+x(9)+x(10);
y4=x(1)+x(2)+x(3)+c(4)*x(4)+x(5)+x(6)+x(7)+x(8)+x(9)+x(10);
y5=x(1)+x(2)+x(3)+x(4)+c(5)*x(5)+x(6)+x(7)+x(8)+x(9)+x(10);
y6=x(1)+x(2)+x(3)+x(4)+x(5)+c(6)*x(6)+x(7)+x(8)+x(9)+x(10);
y7=x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+c(7)*x(7)+x(8)+x(9)+x(10);
y8=x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+c(8)*x(8)+x(9)+x(10);
y9=x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8)+c(9)*x(9)+x(10);
% y10=x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8)+x(9)+c(10)*x(10);

val=[y1,y2,y3,y4,y5,y6,y7,y8 y9]'; 
if ~optimal_policy
    val=[val;y10];
end

%{
V = multivariate_taylor_approximation(@test_system_fun,rand(10,1),2)
%}