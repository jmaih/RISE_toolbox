%--- help for rprt/description ---
%
%  DESCRIPTION Adds a description list to the report.
% 
%  Usage:
%    obj.description(Items)
% 
%  Inputs:
%    - obj: The report object.
%    - Items: A cell array of structs representing the items in the list.
%      Each struct should have 'term' and 'description' fields.
% 
%  Options:
%    None
% 
%  Example:
% <<
%  d=struct('term',{},'description',{});
%  d(1).term='DSGE';
%  d(1).description='Dynamic Stochastic General Equilibrium';
%  d(2).term='RBC';
%  d(2).description='Real Business Cycle';
%  d(3).term='RISE';
%  d(3).description='Rationality In Switching Environments';
%    rprt.description(d);
% >>
%