function B=shuffle_tensor1(A,matsizes,newOrder)

[ra,ca]=size(A);

[~,~,newira,newica]=utils.kronecker.find_reordering([ra,ca],matsizes,newOrder);

B=A(newira,newica);

end