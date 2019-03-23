function B=set_vectorized_mapping(A)
% returns a 3-column matrix in which:
% - The first column represents the rows
% - The second column represents the columns
% - The third column represents the values expander, such that vals(col3)
% returns a vector of same length as the rows of matrix B
%
% A is a cell array of cell arrays
% - in each entry of A is a cell array
%   - where the first element is the row of the matrix
%   - the remaining entries are scalars or vectors that indicate the degree
%   of multiplicity of the corresponding derivative and the location of
%   those derivatives in the final matrix

B=A;

offset=0;

for irow=1:size(A,1)
    
    r=A{irow}{1}; % row
    
    c0=A{irow}(2:end); % vectors of columns
    
    c0=c0(:);
    
    nc=numel(c0);
    
    c=c0;
    
    nr=0;
    
    for ic=1:nc
        
        nitems=numel(c{ic});
        
        c{ic}=ic*ones(nitems,1)+offset;
        
        nr=nr+nitems;
        
    end
    
    offset=c{nc}(end);
    
    c=cell2mat(c);
    
    nr=numel(c);
    
    B{irow}=[r*ones(nr,1),cell2mat(c0),c(:)];
    
end

B=cell2mat(B);

end
