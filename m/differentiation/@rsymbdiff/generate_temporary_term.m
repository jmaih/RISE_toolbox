%--- help for rsymbdiff.generate_temporary_term ---
%
%  GENERATE_TEMPORARY_TERM - Generate a non-conflicting temporary term
% 
%    temporary_term = generate_temporary_term(TMP, existing_symbols) generates
%    a temporary term based on the input TMP (temporary term base) and ensures
%    that the generated term does not conflict with any of the existing_symbols.
% 
%  Input:
%    - TMP (char): The base temporary term.
%    - existing_symbols (cellstr): A cell array of existing symbols to check
%      for conflicts.
% 
%  Output:
%    - temporary_term (char): A temporary term that does not conflict with
%      existing_symbols.
% 
%  Example:
%    TMP = 'temp';
%    existing_symbols = {'x', 'y', 'temp', 'temp1'};
%    temporary_term = generate_temporary_term(TMP, existing_symbols);
% 
%  See also: ismember, num2str
%