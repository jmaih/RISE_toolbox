%--- help for rsymbdiff.writeToFile ---
%
% 
%    WRITE_TO_FILE(d, iter_mirror, inputList, filename) writes symbolic
%    derivatives represented in the 'd' structure as a function to an M-file
%    named 'filename'. The generated function can be used for efficient
%    evaluation of the derivatives.
% 
%    Inputs:
%    - d: Structure containing symbolic derivatives for different equations
%         and orders of approximation.
%    - iter_mirror: Total number of temporary terms.
%    - inputList: List of input variable names.
%    - filename: Name of the generated M-file (without the ".m" extension).
%    - replaceIfelseif (optional): replace occurrences of if_elseif with
%      if_elseiff, so as to avoid evaluating entries that are not used.
% 
%    Outputs:
%    - retcode (=0): flag indicating there is no problem
% 
%    Example:
%    writeToFile(d, iter_mirror, {'x', 'y', 'z'}, 'my_derivatives');
% 
%    See also: rsymbdiff, rsymbdiff.evaluate, rsymbdiff.anonymize
%
%    Other uses of writeToFile
%
%       splanar/writeToFile
%