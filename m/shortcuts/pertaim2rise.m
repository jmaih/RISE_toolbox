%  Simple converter of perturbation AiM models into the RISE platform 
% 
%  [status] = pertaim2rise(aimFileName, riseFileName, parametervalues,loglinearizevars)
%  [status] = pertaim2rise(aimFileName, [], parametervalues,loglinearizevars)
% 
%    Inputs :
% 
%    - aimFileName [char] : name of the perturbationAIM model file with or
%    without extension
% 
%    - riseFileName [char|empty] : Name of the RISE file. If empty, the name
%    will be the same as that of the aim file
% 
%    - parametervalues [char] : Name of the variable containing parameter
%    values in the model file
% 
%    - loglinearizevars [char] : Name of the variable containing the list of
%    the variables to log-linearize in the model file
%