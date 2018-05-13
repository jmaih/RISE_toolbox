function c=find_combos(acell,bcell,c)
% find_combos -- find unique permutations arising from tensors of sums
%
% ::
%
%
%   c=find_combos(acell,bcell)
%
% Args:
%
%    - **acell** [2x1 or 1x2 cell]: where the first element is a number or a
%    char and the second element is an integer describing the number of times
%    the element in the first cell occurs
%
%    - **bcell** [2x1 or 1x2 cell]: where the first element is a number or a
%    char and the second element is an integer describing the number of times
%    the element in the first cell occurs
%
% Returns:
%    :
%
%    - **c** [n x k matrix]: char or doubles depending on the type of inputs
%
% Note:
%
%    - using the Pascal's triangle, we know that
%    (a+b)^5=a^5+5a^4*b+10a^3*b^2+10a^2*b^3+5a*b^4+b^5. Suppose now that a and
%    b are matrices and let's focus on one term, say 10a^3*b^2. We know that
%    we will have 10 combinations of kronecker products in which a appears 3
%    times and b appears twice. But then how to find them? This routine is
%    designed to solve that problem efficiently.
%
% Example:
%
%    - suppose we want to find all permutations in which u appears 4 times and
%    h appears twice. c=find_combos({'u',4},{'h',2}),size(c)
%
%    - suppose we want to find all permutations in which 1 appears 4 times and
%    2 appears three times. c=find_combos({1,4},{2,3})
%
%    See also:
if nargin<3
    
    c=[];

end

num_type=isa(acell{1},'double');

if num_type

    is_error=~isa(bcell{1},'double');

else
    
    is_error=~(ischar(acell{1}) && ischar(bcell{1}));

end

if is_error

    error('the first element in both cell must be a number or a char')

end

cold=c;

nrows=max(size(cold,1),1);

is_update=false;

do_one(acell,bcell)

do_one(bcell,acell)

    function do_one(main,altern)

        if main{2}
        
            co=[cold,main{1}(ones(nrows,1))];
            
            main{2}=main{2}-1;
            
            co=utils.gridfuncs.find_combos(main,altern,co);
            
            if is_update
            
                c=[c
                    co];
            
            else
                
                c=co;
                
                is_update=true;
            
            end
            
        end
        
    end

end

%{
function c=find_combos(a,b,na,nb,c)
if nargin<5
    c=[];
end

rc=max(size(c,1),1);

ca=[];
if na
    co=augment_c(c,a);
    ca=find_combos(a,b,na-1,nb,co);
end
cb=[];
if nb
    co=augment_c(c,b);
    cb=find_combos(a,b,na,nb-1,co);
end

if na||nb
    c=[
        ca
        cb
        ];
end

    function co=augment_c(c,a)
        co=[c,a(ones(rc,1))];
    end

end
%}