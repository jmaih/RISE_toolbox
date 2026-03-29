%--- help for rdico.update_names ---
%
%  update_input_and_output_names modifies names until none of them matches
%  any tokens in the text. One area of application is when for instance we
%  want to create a function for the steady state model and there are
%  variables and/or parameter names that match the input names for the
%  function.
% 
%  Example
% 
%  text = 'This is a sample text, with some names. This text contains some
%  repeated names like John and Jane!'; 
% 
%  List of names to modify
%  names = {'John', 'Jane', 'Bob', 'Alice'};
% 
%  names=rdico.update_names(names,text)
%