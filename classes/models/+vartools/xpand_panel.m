function d=xpand_panel(d)
% going from nvar x T x npages to nvar*npages x T
% where npages is the number of panel members
% the resulting ordering is
% [v1_1,v1_2,...,v1_npages,v2_1,...,v2_npages,...,vnvar_1,...,vnvar_npages]

d=permute(d,[2,3,1]);

d=d(:,:)';

end
