function B=form_B(x,pairs,bsize)
B=zeros(bsize);
B(pairs(:,1))=x(pairs(:,2));
end
