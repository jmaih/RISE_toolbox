function B=shuffle_tensor2(A,matsizes,newOrder)

[ra,ca]=size(A);

[lm,rm]=utils.kronecker.find_reordering([ra,ca],matsizes,newOrder);

B=lm*A*rm;

end