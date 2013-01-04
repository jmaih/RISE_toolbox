function [H,this,hjmap,jac_hess_mat]=hessian(funcs,varList,wrt)

if ischar(funcs)
    funcs=cellstr(funcs);
elseif iscell(funcs)
    if isa(funcs{1},'sad_reverse')
        funcs=[funcs{:}];
    end
end
if ischar(wrt)
    wrt=cellstr(wrt);
end
eqtn_nbr=numel(funcs);
wrt_nbr=numel(wrt);

if eqtn_nbr>1
    error('number of equations must be 1')
end

% call the jacobian in non-vectorized form
[this,jmap,jac_restrict_refs]=sad_reverse.jacobian(funcs,varList,wrt,false);

% initialize the Hessian matrix
H=sad_reverse(0);
H=H(ones(wrt_nbr));

for irow=1:wrt_nbr
    % call the rest in non-vectorized form as well
    H(irow:end,irow)=sad_reverse.jacobian(this(irow),varList,wrt(irow:end),false);
end

jmap.fixed_references=union(jmap.fixed_references,...
    sad_reverse.extract_list(jac_restrict_refs,jmap.prefix_list(1)));
% build on the existing map to create something for the hessian
[hess_restrict_refs,hmap]=char(H,jmap);
hjmap=hmap;
% cut the head as it belongs to jmap
hmap.fid=hmap.fid(jmap.line_count+1:end,:);
hmap.line_count=size(hmap.fid,1);
hjmap.fixed_references=union(jmap.fixed_references,...
    sad_reverse.extract_list(hess_restrict_refs,hmap.prefix_list(1)));
hmat=sad_reverse.create_matrix(hmap,hess_restrict_refs,wrt_nbr,'Hess_');
% closing statement to include the elements above the diagonal
hmat=[hmat
    'Hess_=tril(Hess_)+tril(Hess_,-1)'';'];
% create a matrix for the jacobian as well
jmat=sad_reverse.create_matrix(jmap,jac_restrict_refs,wrt_nbr,'Jac_');

% put jacobian and hessian matrices together
jac_hess_mat=[jmat;hmat];

