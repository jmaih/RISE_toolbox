function P=var_rotation(fa0ap,Q)
n=size(fa0ap,2);
[Q,~,tags]=sort_Q(Q);
P=nan(n);
for jj=1:n
    Qj=Q{1,tags{jj}};
    Qtilde=[Qj*fa0ap
        P(:,1:jj-1)'];
    [q,~]=qr(Qtilde');
    P(:,jj)=q(end,:)';
end
% P(:,tags)=P;