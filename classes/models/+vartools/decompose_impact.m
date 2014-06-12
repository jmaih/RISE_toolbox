function [omg,sig]=decompose_impact(C)

big=max(abs(C),[],1);
sig=diag(big);
omg=C./repmat(big,size(C,1),1);