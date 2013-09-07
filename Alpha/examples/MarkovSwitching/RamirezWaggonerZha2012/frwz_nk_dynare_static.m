function [residual, g1, g2] = frwz_nk_dynare_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 3, 1);

%
% Model equations
%

T16 = 1-.5*params(3)*(y(1)-1)^2;
T48 = (y(3)/(y(3)))^params(7);
T54 = (y(1)/(y(1)))^((1-params(7))*params(6));
residual(1) = 1-params(1)*T16*y(2)*y(3)/(y(1)*T16*y(2)*exp(params(4)));
residual(2) = 1-params(2)+y(2)*T16*params(2)+y(1)*(y(1)-1)*T16*params(1)*params(3)/T16-y(1)*params(3)*(y(1)-1);
residual(3) = T48*T54*exp(params(8)*x(1))-y(3)/(y(3));
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(3, 3);

  %
  % Jacobian matrix
  %

  g1(1,1)=(-((y(1)*T16*y(2)*exp(params(4))*y(3)*y(2)*params(1)*(-(.5*params(3)*2*(y(1)-1)))-params(1)*T16*y(2)*y(3)*(T16*y(2)*exp(params(4))+y(1)*exp(params(4))*y(2)*(-(.5*params(3)*2*(y(1)-1)))))/(y(1)*T16*y(2)*exp(params(4))*y(1)*T16*y(2)*exp(params(4)))));
  g1(1,2)=(-((y(1)*T16*y(2)*exp(params(4))*params(1)*T16*y(3)-params(1)*T16*y(2)*y(3)*y(1)*T16*exp(params(4)))/(y(1)*T16*y(2)*exp(params(4))*y(1)*T16*y(2)*exp(params(4)))));
  g1(1,3)=(-(params(1)*T16*y(2)/(y(1)*T16*y(2)*exp(params(4)))));
  g1(2,1)=y(2)*params(2)*(-(.5*params(3)*2*(y(1)-1)))+(T16*((y(1)-1)*T16*params(1)*params(3)+y(1)*(T16*params(1)*params(3)+(y(1)-1)*params(1)*params(3)*(-(.5*params(3)*2*(y(1)-1)))))-y(1)*(y(1)-1)*T16*params(1)*params(3)*(-(.5*params(3)*2*(y(1)-1))))/(T16*T16)-(params(3)*(y(1)-1)+params(3)*y(1));
  g1(2,2)=T16*params(2);
  g1(3,1)=exp(params(8)*x(1))*T48*((y(1))-y(1))/((y(1))*(y(1)))*getPowerDeriv(y(1)/(y(1)),(1-params(7))*params(6),1);
  g1(3,3)=exp(params(8)*x(1))*T54*((y(3))-y(3))/((y(3))*(y(3)))*getPowerDeriv(y(3)/(y(3)),params(7),1)-((y(3))-y(3))/((y(3))*(y(3)));
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],3,9);
end
end
