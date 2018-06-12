function [lia,locb] = ismember(A,B)
% INTERNAL FUNCTION: checks membership of the rows of a matrix in another matrix
%
% ::
%
%   [lia] = ismember(A,B)
%   [lia,locb] = ismember(A,B)
%
% Args:
%
%    - **A** [n x k matrix] : matrix whose rows we want to locate or check the
%      membership
%    - **B** [m x k matrix] : reference matrix
%
% Returns:
%    :
%
%    - **lia** [n x 1 logical] : true for rows that belong
%    - **locb** [n x 1 vector] : location of the elements of A in B
%

if size(A,2)~=size(B,2)

    error('matrices should have the same number of columns')

end

% Duplicates within the sets are eliminated
[uA,~,icA] = unique(A,'rows','sorted');

[uB,ib] = unique(B,'rows','sorted');

% Sort the unique elements of A and B, duplicate entries are adjacent
[sortuAuB,IndSortuAuB] = sortrows([uA;uB]);

% Find matching entries
d = sortuAuB(1:end-1,:)==sortuAuB(2:end,:);     % d indicates matching entries

d = all(d,2);                                   % Finds the index of matching entries

ndx1 = IndSortuAuB(d);                          % NDX1 are locations of repeats in C

szuA = size(uA,1);
% Find locb by using given indices
[ndx1,idx] = sort(ndx1(:));

[lia, locb] = builtin('_ismemberhelper',icA,ndx1);

locb(lia) = idx(locb(lia));

d = find(d);

newd = d(locb(lia));                    % NEWD is D for non-unique A

where = ib(IndSortuAuB(newd+1)-szuA);   % Index values of uB through UNIQUE

locb(lia) = where;                      % Return first or last occurrence of A within B

end
