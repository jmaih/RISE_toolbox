function M=set_up_moments(order,nshocks,mom_vector)

% mom_vector can be
% - a function handle: it returns a vector containing the moments at each
% order or a matrix, in which case the columns represent the orders and the
% rows the different shocks
% - a vector: the number of elements corresponds to the order of the
% moments.
% - a matrix, in which case the columns represent the orders and the
% rows the different shocks 
% mom_vector is a vector or produces a vector, it is assumed that all
% shocks have the same distribution

% hops.set_up_moments(5,7,cell2mat(hops.set_up_normal_moments(5,1)))

% hops.set_up_moments(5,7,@hops.set_up_normal_moments)

if isa(mom_vector,'function_handle')
    
    mom_vector=cell2mat(mom_vector(order,1));
    
end

if isvector(mom_vector)
    
    mom_vector=mom_vector(:).';
    
    mom_vector=mom_vector(ones(1,nshocks),:);
    
end

if ~isequal(size(mom_vector),[nshocks,order])
    
    error('wrong size of mom_vector')
    
end

M=cell(1,order);

for oo=1:order
    
    M{oo}=set_up_one(oo);
    
end

    function m=set_up_one(o)
        
        these_moms=mom_vector(:,o);
        
        if o==1
            
            m=these_moms;
            
        else
            
            siz=nshocks*ones(1,o);
            
            m=zeros(siz);
            
            e=repmat({(1:nshocks).'},1,o);
            
            IND=sub2ind(siz,e{:});
            
            m(IND)=these_moms;
            
        end
        
    end

end
