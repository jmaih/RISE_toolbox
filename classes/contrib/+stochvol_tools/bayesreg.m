function [bdraw,b1,S1] = bayesreg(b0,S0,sv,y,x)

% bdraw = bayesreg[b0,S0,sv,y,x];
% bayesian regression
% b0 = prior mean for regression parameters
% S0 = prior variance for regression parameters
% sv = standard error for regression errors, assumed known
% y = Tx1 vector of dependent variables
% x = Txk matrix of independent variables
% bdraw is a draw from the posterior

S1 = inv(inv(S0) + x'*x/(sv^2));
b1 = S1*(inv(S0)*b0 + x'*y/(sv^2));
%bdraw = b1 + sqrtm(S1)*randn(size(b1));
%bdraw = b1 + cholfortran(S1)*randn(size(b1));
bdraw = b1 + chol(S1)*randn(size(b1));