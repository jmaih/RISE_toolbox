%--- help for rdico.cell2matize ---
%
%  CELL2MATIZE - Convert a cell array of strings into a formatted string
% 
%    c2m = cell2matize(list) converts a cell array of strings 'list' into a
%    formatted string where the strings are joined using the default separator
%    '|', and enclosed within default brackets '()'.
% 
%    c2m = cell2matize(list, separator) allows you to specify a custom
%    separator character to join the strings in the list.
% 
%    c2m = cell2matize(list, separator, brackets) allows you to specify custom
%    opening and closing brackets for the formatted string.
% 
%  Input:
%    - list (cellstr): A cell array of strings to be formatted.
%    - separator (char | {'|'}): (Optional) The character used to separate the
%      strings in the formatted output.
%    - brackets (cell of two chars | {'()'}): (Optional) A cell array of two
%      characters specifying the opening and closing brackets for the formatted
%      output.
% 
%  Output:
%    - c2m (char): The formatted string obtained by joining the input strings
%      with the specified separator and enclosing them within brackets.
% 
%  Examples:
%    c2m = cell2matize({'x','y','z','t'});
%    c2m = cell2matize({'x','y','z','t'}, ',');
%    c2m = cell2matize({'x','y','z','t'}, [], {'[',']'});
% 
%  See also: strcat, cell2mat
%