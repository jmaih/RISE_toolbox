function [lm,rm,newira,newica]=find_reordering(siza,matsizes,newOrder)

% siza : size of the entire kronecker product

% matsizes : sizes of the matrices entering the "ideal" kronecker product

% newOrder: desired order of the matrices whose sizes are in matsizes

% lm, rm sparse matrices such that B=lm*A*rm gives the kronecker product
% with desired reshuffling.

ra=siza(1); ca=siza(2);

rows=matsizes(end:-1:1,1);

cols=matsizes(end:-1:1,2);

ira=reshape(1:ra,rows(:).');

ica=reshape(1:ca,cols(:).');

nmat=numel(rows);

order_=fliplr(1:nmat);

newOrder_=fliplr(newOrder(:).');

reorder=nan(1,nmat);

for ii=1:nmat
    
    reorder(ii)=find(order_==newOrder_(ii));
    
end

newira=permute(ira,reorder);

newica=permute(ica,reorder);

newira=newira(:);

newica=newica(:);

posl(newira)=1:ra;

posr(newica)=1:ca;

lm=speye(ra); lm=lm(:,posl);

rm=speye(ca); rm=rm(posr,:);

end
