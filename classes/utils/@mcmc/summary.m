%--- help for categorical/summary ---
%
% SUMMARY Print summary of a categorical array.
%    SUMMARY(A) displays the number of elements in the categorical array A
%    that are equal to each of A's categories.  If A contains any undefined
%    elements, the output also includes the number of undefined elements.
% 
%    If A is a vector, SUMMARY displays counts for the entire vector.  If A is a
%    matrix or N-D array, SUMMARY displays counts separately for each column of A.
% 
%    SUMMARY(A,DIM) displays the summary computed along the dimension DIM of A.
% 
%    See also ISCATEGORY, ISMEMBER, COUNTCATS.
%
%    Documentation for categorical/summary
%       doc categorical/summary
%
%    Other uses of summary
%
%       backtestEngine/summary               distributed/summary
%       brinsonAttribution/summary           DriftDiagnostics/summary
%       clibgen.LibraryDefinition/summary    mcmc/summary
%       codistributed/summary                tabular/summary
%       dataset/summary                      tall/summary
%