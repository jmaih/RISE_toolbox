%--- help for struct/merge ---
%
%  Merge structures with unique fields.
% 
%  Syntax:
%    A = merge(struct1, struct2, ...)
% 
%  Description:
%    The `merge` function takes multiple structures as input and merges
%    them into a single structure, ensuring that the resulting structure
%    has unique field names. It concatenates the data from each input
%    structure.
% 
%  Input:
%    - varargin: Input structures to be merged.
% 
%  Output:
%    - A: Merged structure.
% 
%  Example:
%    s1.field1 = 1;
%    s1.field2 = 'hello';
%    s2.field2 = 'hello';
%    s2.field3 = [1 2 3];
%    mergedStruct = merge(s1, s2);
%    disp(mergedStruct);
% 
%    % Output:
%    %   field1: 1
%    %   field2: 'hello'
%    %   field3: [1 2 3]
% 
%  Notes:
%    - If field names are not unique, they should have the same values in
%      order to avoid errors.
% 
%  See also:
%    - struct, orderfields
%
%    Documentation for struct/merge
%       doc struct/merge
%
%