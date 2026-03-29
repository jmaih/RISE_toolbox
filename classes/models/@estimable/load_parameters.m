%--- help for estimable/load_parameters ---
%
%  load_parameters loads the parameters. This allows the user to quickly
%  load the parameters from a file, which may be the output of estimation, and get
%  going with irfs, simulations, etc.
% 
%  ::
% 
%     model=load_parameters(model,the_mode_file)
% 
%  Args:
% 
%     - **model** (estimable object): model object
% 
%     - **the_mode_file** (m-file): file containing the parameters and their
%       values
% 
%  Returns:
% 
%     - **model** (estimable object): reparameterized model object
%