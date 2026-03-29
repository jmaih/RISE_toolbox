%--- help for struct/renamefields ---
%
%  RENAMESFIELDS Renames fields in a structure
% 
%    newStruct = renamefields(oldStruct, oldNames, newNames) returns a
%    structure with the same data as oldStruct, but with fields named
%    oldNames renamed to newNames.
% 
%    Example:
%        oldStruct = struct('foo', 1, 'bar', 2);
%        newStruct = renamefields(oldStruct, {'foo', 'bar'}, {'alpha', 'beta'});
%        % newStruct is struct('alpha', 1, 'beta', 2)
%
%    Documentation for struct/renamefields
%       doc struct/renamefields
%
%