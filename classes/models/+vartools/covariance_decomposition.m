function [D,W]=covariance_decomposition(SIGu)
% This function decomposes SIGu = W*diag(D)^2*W', where D is a
% vector of positive numbers and W is lower triangular matrix with
% ones on its main diagonal.
% see e.g. Lutkepohl(2005):"New Introduction to Multiple
% Time-Series Analysis" page 58
P=chol(SIGu,'lower');
D=diag(P);
W=bsxfun(@times,P,1./D.');
end
