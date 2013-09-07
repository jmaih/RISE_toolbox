function [residual, g1, g2] = constant_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 13, 1);

%
% Model equations
%

paibar__ = params(9);
T40 = y(2)^(1-params(6));
T52 = (1-params(6))^(1-params(6));
T55 = y(1)^(1-params(6))*y(10)^params(6)*T52*params(6)^params(6);
T66 = y(6)/paibar__-1;
T69 = y(6)*params(7)*T66/paibar__;
T99 = (y(6)/(y(6)))^params(10);
T103 = (y(8)/(y(8)))^params(11);
T200 = 1/params(12)/(y(7)/params(12));
lhs =y(1);
rhs =params(1)*y(2)^params(2)*y(3)^params(3);
residual(1)= lhs-rhs;
lhs =1;
rhs =y(4)*y(5)/y(6);
residual(2)= lhs-rhs;
lhs =1;
rhs =params(4)*(1+y(10)-params(5));
residual(3)= lhs-rhs;
lhs =y(8);
rhs =y(7)*y(12)^params(6)*T40;
residual(4)= lhs-rhs;
lhs =y(12);
rhs =y(12)*(1-params(5))+y(11);
residual(5)= lhs-rhs;
lhs =y(9);
rhs =T55/y(7);
residual(6)= lhs-rhs;
lhs =y(12)/y(2);
rhs =y(1)*params(6)/(1-params(6))/y(10);
residual(7)= lhs-rhs;
lhs =T69;
rhs =1-params(8)+y(9)*params(8)+y(8)*y(6)*T66*y(5)*params(7)/paibar__/y(8);
residual(8)= lhs-rhs;
lhs =y(3)+y(11)+params(18)*(y(8));
rhs =y(8)*(1-params(7)*0.5*T66^2);
residual(9)= lhs-rhs;
lhs =y(13);
rhs =(y(4))*T99*T103;
residual(10)= lhs-rhs;
lhs =y(4);
rhs =y(13);
residual(11)= lhs-rhs;
lhs =log(y(7)/params(12));
rhs =log(y(7)/params(12))*params(13)+params(16)*x(1);
residual(12)= lhs-rhs;
lhs =y(5);
rhs =params(4)*(y(5)/params(4))^params(14)*exp(params(15)*x(2));
residual(13)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(13, 13);

  %
  % Jacobian matrix
  %

  g1(1,1)=1;
  g1(1,2)=(-(y(3)^params(3)*params(1)*getPowerDeriv(y(2),params(2),1)));
  g1(1,3)=(-(params(1)*y(2)^params(2)*getPowerDeriv(y(3),params(3),1)));
  g1(2,4)=(-(y(5)/y(6)));
  g1(2,5)=(-(y(4)/y(6)));
  g1(2,6)=(-((-(y(4)*y(5)))/(y(6)*y(6))));
  g1(3,10)=(-params(4));
  g1(4,2)=(-(y(7)*y(12)^params(6)*getPowerDeriv(y(2),1-params(6),1)));
  g1(4,7)=(-(y(12)^params(6)*T40));
  g1(4,8)=1;
  g1(4,12)=(-(T40*y(7)*getPowerDeriv(y(12),params(6),1)));
  g1(5,11)=(-1);
  g1(5,12)=1-(1-params(5));
  g1(6,1)=(-(params(6)^params(6)*T52*y(10)^params(6)*getPowerDeriv(y(1),1-params(6),1)/y(7)));
  g1(6,7)=(-((-T55)/(y(7)*y(7))));
  g1(6,9)=1;
  g1(6,10)=(-(params(6)^params(6)*T52*y(1)^(1-params(6))*getPowerDeriv(y(10),params(6),1)/y(7)));
  g1(7,1)=(-(params(6)/(1-params(6))/y(10)));
  g1(7,2)=(-y(12))/(y(2)*y(2));
  g1(7,10)=(-((-(y(1)*params(6)/(1-params(6))))/(y(10)*y(10))));
  g1(7,12)=1/y(2);
  g1(8,5)=(-(y(8)*T69/y(8)));
  g1(8,6)=(params(7)*T66+y(6)*params(7)*1/paibar__)/paibar__-y(8)*(T66*y(5)*params(7)+y(6)*y(5)*params(7)*1/paibar__)/paibar__/y(8);
  g1(8,9)=(-params(8));
  g1(9,3)=1;
  g1(9,6)=(-(y(8)*(-(params(7)*0.5*1/paibar__*2*T66))));
  g1(9,8)=params(18)-(1-params(7)*0.5*T66^2);
  g1(9,11)=1;
  g1(10,4)=(-(T99*T103));
  g1(10,6)=(-(T103*(y(4))*((y(6))-y(6))/((y(6))*(y(6)))*getPowerDeriv(y(6)/(y(6)),params(10),1)));
  g1(10,8)=(-((y(4))*T99*((y(8))-y(8))/((y(8))*(y(8)))*getPowerDeriv(y(8)/(y(8)),params(11),1)));
  g1(10,13)=1;
  g1(11,4)=1;
  g1(11,13)=(-1);
  g1(12,7)=T200-params(13)*T200;
  g1(13,5)=1-exp(params(15)*x(2))*params(4)*1/params(4)*getPowerDeriv(y(5)/params(4),params(14),1);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],13,169);
end
end
