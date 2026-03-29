% ~=  Not equal.
%    A ~= B does element by element comparisons between A and B and returns
%    an array with elements set to logical 1 (TRUE) where the relation is
%    true and elements set to logical 0 (FALSE) where it is not. A and B
%    must have compatible sizes. In the simplest cases, they can be the same
%    size or one can be a scalar. Two inputs have compatible sizes if, for
%    every dimension, the dimension sizes of the inputs are either the same
%    or one of them is 1.
% 
%    C = NE(A,B) is called for the syntax 'A ~= B' when A or B is an object.
% 
%    See <a href="matlab:helpview('matlab','MATLAB_OPS')">MATLAB Operators and Special Characters</a> for more details.
%
%    Documentation for ne
%       doc ne
%
%    Other uses of ne
%
%       adolm/ne            mtree/ne               splanar/ne
%       categorical/ne      opaque/ne              string/ne
%       codistributed/ne    qrandstream/ne         sym/ne
%       datetime/ne         rise_dates.dates/ne    symbolic/ne
%       duration/ne         rsymbdiff/ne           tabular/ne
%       gpuArray/ne         scribehandle/ne        timer/ne
%       handle/ne           serial/ne              ts/ne
%       instrument/ne
%