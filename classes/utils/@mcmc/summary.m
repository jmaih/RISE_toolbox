%--- help for tabular/summary ---
%
% SUMMARY Print summary of a table or a timetable.
%    SUMMARY(T) prints a summary of T and the variables that it contains. If
%    T is a table, then SUMMARY displays the description from
%    T.Properties.Description followed by a summary of the table variables.
%    If T is a timetable, then SUMMARY additionally displays a summary of
%    the row times.
% 
%    S = SUMMARY(T) returns a summary of T as a structure. Each field of S
%    contains a summary for the corresponding variable in T. If T is a
%    timetable, S contains an additional field for a summary of the row
%    times.
%
%    Other functions named summary
%
%       categorical/summary      dataset/summary        mcmc/summary
%       codistributed/summary    distributed/summary    tall/summary
%