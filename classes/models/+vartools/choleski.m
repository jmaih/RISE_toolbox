function [Rfunc,ident]=choleski(n,ordering)

if nargin < 2
    
    ordering=[];
    
end

if isempty(ordering)
    
    ordering=struct('old',1:n,'new',1:n);
    
end

old_order=ordering.old;

new_order=ordering.new;

if ischar(old_order),old_order=cellstr(old_order); end

if ischar(new_order),new_order=cellstr(new_order); end

if iscellstr(old_order) && ~iscellstr(new_order)
    
    error('both new_order and old_order should be of the same type')
    
end

nitems=numel(old_order);

if ~(all(ismember(new_order,old_order)) &&...
        numel(new_order)==nitems)
    
    error('both new_order and old_order should have the same elements')
    
end

if nitems ~=n
    
    error('number of variables does not match the size of the system')
    
end

order=1:n;

if iscellstr(new_order)
    
    for iv=1:n
        
        order(iv)=find(strcmp(new_order{iv},old_order));
        
    end
    
elseif isnumeric(new_order) && ...
        all(ismember(new_order,order))
    
    order=new_order;
    
else
    
    error('wrong input format')
    
end

iorder(order)=1:n;

ident = struct('isexact',true,'isover',false,'isunder',false);

Rfunc=@choleski_engine;

    function [R,retcode]=choleski_engine(p)
         
        retcode=0;
        
        R = chol(p.S(order,order),'lower');
        
        R = R(iorder,iorder);
        
    end


end