function [equation_blocks,variable_blocks]=block_triangularize(Incidence,do_blocks)
% block_triangularize -- creates blocks of equation-variables so as to
% permit a recursive solution
%
% ::
%
%
%   [equation_blocks,variable_blocks]=block_triangularize(Incidence,do_blocks)
%
% Args:
%
%    - **Incidence** [sparse]: Incidence matrix where rows represent equations
%    and columns represent variables.
%
%    - **do_blocks** [false|{true}]: triggers the computation of blocks.
%
% Returns:
%    :
%
%    - **equation_blocks** [cell array]: blocks of equations IDs
%
%    - **variable_blocks** [cell array]: blocks of variables IDs
%
% Note:
%
% Example:
%
%    See also:

% References
%-------------
% Pothen, Alex and Chin-Ju Fan "Computing the Block Triangular Form of a
% Sparse Matrix" ACM Transactions on Mathematical Software Vol 16, No. 4
% Dec. 1990, pp. 303-324.  
%
% Inspired by Jaromir Benes' Blazer

if nargin<2
    do_blocks=true;
end

Incidence=sparse(Incidence);

[r,c]=size(Incidence);
equations=1:r;
variables=1:c;

if do_blocks
    [ordered_incidence, ordered_equations, ordered_variables] = reorder();
    %  this.IsSingular = sprank(ordered_incidence)<min(size(ordered_incidence));
    [equation_blocks,variable_blocks] = make_blocks();
else
    equation_blocks = {equations};
    variable_blocks = {variables};
end

    function [equation_blocks,variable_blocks] = make_blocks()
        equation_blocks = cell(1,0);
        variable_blocks = cell(1,0);
        currEqtnBlk = [];
        currNameBlk = [];
        for irow = r : -1 : 1
            currNameBlk(end+1) = ordered_variables(irow); %#ok<AGROW>
            currEqtnBlk(end+1) = ordered_equations(irow); %#ok<AGROW>
            if ~any(any(ordered_incidence(irow:end, 1:irow-1)))
                variable_blocks{end+1} = currNameBlk; %#ok<AGROW>
                equation_blocks{end+1} = currEqtnBlk; %#ok<AGROW>
                currNameBlk = [];
                currEqtnBlk = [];
            end
        end
    end

    function [ordered_incidence,ordered_equations,ordered_variables] = reorder()
        order_cols = colamd(Incidence);
        ordered_incidence = Incidence(:,order_cols);
        ordered_variables = variables(order_cols);
        
        order_rows = colamd(ordered_incidence.');
        ordered_incidence = ordered_incidence(order_rows,:);
        ordered_equations = equations(order_rows);
        
        [order_rows2,order_cols2] = dmperm(ordered_incidence);
        ordered_incidence = ordered_incidence(order_rows2,order_cols2);
        ordered_variables = ordered_variables(order_cols2);
        ordered_equations = ordered_equations(order_rows2);
    end

end
