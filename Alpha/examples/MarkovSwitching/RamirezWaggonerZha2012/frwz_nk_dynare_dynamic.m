function [residual, g1, g2, g3] = frwz_nk_dynare_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(3, 1);
T16 = 1-.5*params(3)*(y(2)-1)^2;
T26 = 1-.5*params(3)*(y(5)-1)^2;
T31 = T26*y(6)*exp(params(4));
T32 = y(5)*T31;
T55 = (y(1)/(steady_state(3)))^params(7);
T61 = (y(2)/(steady_state(1)))^((1-params(7))*params(6));
T91 = 1/(steady_state(1))*getPowerDeriv(y(2)/(steady_state(1)),(1-params(7))*params(6),1);
T100 = T31+y(5)*exp(params(4))*y(6)*(-(.5*params(3)*2*(y(5)-1)));
T110 = T26*((y(5)-1)*T16*params(1)*params(3)+y(5)*T16*params(1)*params(3))-y(5)*(y(5)-1)*T16*params(1)*params(3)*(-(.5*params(3)*2*(y(5)-1)));
T122 = 1/(steady_state(3));
T124 = T122*getPowerDeriv(y(1)/(steady_state(3)),params(7),1);
residual(1) = 1-params(1)*T16*y(3)*y(4)/T32;
residual(2) = 1-params(2)+y(3)*T16*params(2)+y(5)*(y(5)-1)*T16*params(1)*params(3)/T26-y(2)*params(3)*(y(2)-1);
residual(3) = T55*T61*exp(params(8)*x(it_, 1))-y(4)/(steady_state(3));
if nargout >= 2,
  g1 = zeros(3, 7);

  %
  % Jacobian matrix
  %

  g1(1,2)=(-(y(4)*y(3)*params(1)*(-(.5*params(3)*2*(y(2)-1)))/T32));
  g1(1,5)=(-((-(params(1)*T16*y(3)*y(4)*T100))/(T32*T32)));
  g1(1,3)=(-(params(1)*T16*y(4)/T32));
  g1(1,6)=(-((-(params(1)*T16*y(3)*y(4)*y(5)*T26*exp(params(4))))/(T32*T32)));
  g1(1,4)=(-(params(1)*T16*y(3)/T32));
  g1(2,2)=y(3)*params(2)*(-(.5*params(3)*2*(y(2)-1)))+y(5)*(y(5)-1)*params(1)*params(3)*(-(.5*params(3)*2*(y(2)-1)))/T26-(params(3)*(y(2)-1)+params(3)*y(2));
  g1(2,5)=T110/(T26*T26);
  g1(2,3)=T16*params(2);
  g1(3,2)=exp(params(8)*x(it_, 1))*T55*T91;
  g1(3,1)=exp(params(8)*x(it_, 1))*T61*T124;
  g1(3,4)=(-T122);
  g1(3,7)=T55*T61*params(8)*exp(params(8)*x(it_, 1));
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  v2 = zeros(38,3);
  v2(1,1)=1;
  v2(1,2)=9;
  v2(1,3)=(-(y(4)*y(3)*params(1)*(-(2*.5*params(3)))/T32));
  v2(2,1)=1;
  v2(2,2)=30;
  v2(2,3)=(-((-(y(4)*y(3)*params(1)*(-(.5*params(3)*2*(y(2)-1)))*T100))/(T32*T32)));
  v2(3,1)=1;
  v2(3,2)=12;
  v2(3,3)=  v2(2,3);
  v2(4,1)=1;
  v2(4,2)=33;
  v2(4,3)=(-((T32*T32*(-(params(1)*T16*y(3)*y(4)*(exp(params(4))*y(6)*(-(.5*params(3)*2*(y(5)-1)))+exp(params(4))*y(6)*(-(.5*params(3)*2*(y(5)-1)))+y(5)*exp(params(4))*y(6)*(-(2*.5*params(3))))))-(-(params(1)*T16*y(3)*y(4)*T100))*(T32*T100+T32*T100))/(T32*T32*T32*T32)));
  v2(5,1)=1;
  v2(5,2)=16;
  v2(5,3)=(-(y(4)*params(1)*(-(.5*params(3)*2*(y(2)-1)))/T32));
  v2(6,1)=1;
  v2(6,2)=10;
  v2(6,3)=  v2(5,3);
  v2(7,1)=1;
  v2(7,2)=19;
  v2(7,3)=(-((-(T100*params(1)*T16*y(4)))/(T32*T32)));
  v2(8,1)=1;
  v2(8,2)=31;
  v2(8,3)=  v2(7,3);
  v2(9,1)=1;
  v2(9,2)=37;
  v2(9,3)=(-((-(y(4)*y(3)*params(1)*(-(.5*params(3)*2*(y(2)-1)))*y(5)*T26*exp(params(4))))/(T32*T32)));
  v2(10,1)=1;
  v2(10,2)=13;
  v2(10,3)=  v2(9,3);
  v2(11,1)=1;
  v2(11,2)=40;
  v2(11,3)=(-((T32*T32*(-(params(1)*T16*y(3)*y(4)*(T26*exp(params(4))+y(5)*exp(params(4))*(-(.5*params(3)*2*(y(5)-1))))))-(-(params(1)*T16*y(3)*y(4)*y(5)*T26*exp(params(4))))*(T32*T100+T32*T100))/(T32*T32*T32*T32)));
  v2(12,1)=1;
  v2(12,2)=34;
  v2(12,3)=  v2(11,3);
  v2(13,1)=1;
  v2(13,2)=38;
  v2(13,3)=(-((-(params(1)*T16*y(4)*y(5)*T26*exp(params(4))))/(T32*T32)));
  v2(14,1)=1;
  v2(14,2)=20;
  v2(14,3)=  v2(13,3);
  v2(15,1)=1;
  v2(15,2)=41;
  v2(15,3)=(-((-((-(params(1)*T16*y(3)*y(4)*y(5)*T26*exp(params(4))))*(T32*y(5)*T26*exp(params(4))+T32*y(5)*T26*exp(params(4)))))/(T32*T32*T32*T32)));
  v2(16,1)=1;
  v2(16,2)=23;
  v2(16,3)=(-(y(3)*params(1)*(-(.5*params(3)*2*(y(2)-1)))/T32));
  v2(17,1)=1;
  v2(17,2)=11;
  v2(17,3)=  v2(16,3);
  v2(18,1)=1;
  v2(18,2)=26;
  v2(18,3)=(-((-(params(1)*T16*y(3)*T100))/(T32*T32)));
  v2(19,1)=1;
  v2(19,2)=32;
  v2(19,3)=  v2(18,3);
  v2(20,1)=1;
  v2(20,2)=24;
  v2(20,3)=(-(params(1)*T16/T32));
  v2(21,1)=1;
  v2(21,2)=18;
  v2(21,3)=  v2(20,3);
  v2(22,1)=1;
  v2(22,2)=27;
  v2(22,3)=(-((-(params(1)*T16*y(3)*y(5)*T26*exp(params(4))))/(T32*T32)));
  v2(23,1)=1;
  v2(23,2)=39;
  v2(23,3)=  v2(22,3);
  v2(24,1)=2;
  v2(24,2)=9;
  v2(24,3)=y(3)*params(2)*(-(2*.5*params(3)))+y(5)*(y(5)-1)*params(1)*params(3)*(-(2*.5*params(3)))/T26-(params(3)+params(3));
  v2(25,1)=2;
  v2(25,2)=30;
  v2(25,3)=(T26*((y(5)-1)*params(1)*params(3)*(-(.5*params(3)*2*(y(2)-1)))+y(5)*params(1)*params(3)*(-(.5*params(3)*2*(y(2)-1))))-y(5)*(y(5)-1)*params(1)*params(3)*(-(.5*params(3)*2*(y(2)-1)))*(-(.5*params(3)*2*(y(5)-1))))/(T26*T26);
  v2(26,1)=2;
  v2(26,2)=12;
  v2(26,3)=  v2(25,3);
  v2(27,1)=2;
  v2(27,2)=33;
  v2(27,3)=(T26*T26*((-(.5*params(3)*2*(y(5)-1)))*((y(5)-1)*T16*params(1)*params(3)+y(5)*T16*params(1)*params(3))+T26*(T16*params(1)*params(3)+T16*params(1)*params(3))-((-(.5*params(3)*2*(y(5)-1)))*((y(5)-1)*T16*params(1)*params(3)+y(5)*T16*params(1)*params(3))+y(5)*(y(5)-1)*T16*params(1)*params(3)*(-(2*.5*params(3)))))-T110*(T26*(-(.5*params(3)*2*(y(5)-1)))+T26*(-(.5*params(3)*2*(y(5)-1)))))/(T26*T26*T26*T26);
  v2(28,1)=2;
  v2(28,2)=16;
  v2(28,3)=params(2)*(-(.5*params(3)*2*(y(2)-1)));
  v2(29,1)=2;
  v2(29,2)=10;
  v2(29,3)=  v2(28,3);
  v2(30,1)=3;
  v2(30,2)=9;
  v2(30,3)=exp(params(8)*x(it_, 1))*T55*1/(steady_state(1))*1/(steady_state(1))*getPowerDeriv(y(2)/(steady_state(1)),(1-params(7))*params(6),2);
  v2(31,1)=3;
  v2(31,2)=2;
  v2(31,3)=exp(params(8)*x(it_, 1))*T91*T124;
  v2(32,1)=3;
  v2(32,2)=8;
  v2(32,3)=  v2(31,3);
  v2(33,1)=3;
  v2(33,2)=1;
  v2(33,3)=exp(params(8)*x(it_, 1))*T61*T122*T122*getPowerDeriv(y(1)/(steady_state(3)),params(7),2);
  v2(34,1)=3;
  v2(34,2)=44;
  v2(34,3)=T55*T91*params(8)*exp(params(8)*x(it_, 1));
  v2(35,1)=3;
  v2(35,2)=14;
  v2(35,3)=  v2(34,3);
  v2(36,1)=3;
  v2(36,2)=43;
  v2(36,3)=T61*T124*params(8)*exp(params(8)*x(it_, 1));
  v2(37,1)=3;
  v2(37,2)=7;
  v2(37,3)=  v2(36,3);
  v2(38,1)=3;
  v2(38,2)=49;
  v2(38,3)=T55*T61*params(8)*params(8)*exp(params(8)*x(it_, 1));
  g2 = sparse(v2(:,1),v2(:,2),v2(:,3),3,49);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],3,343);
end
end
