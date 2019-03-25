%  reorder_outputs : reorders the outputs of a function
% 
%  ::
% 
%    varargout=reorder_outputs(f,order,varargin)
% 
%  Args:
% 
%     f (function handle): function for which one wants to reorder the
%       outputs
% 
%     order (vector): order in which the outputs are to be reordered. Note
%       that all elements in the reordering should be in 1:numel(order)
% 
%     varargin (varargin): input arguments to the f function
% 
%  Returns:
%     :
% 
%     - **varargout** [varargout]: output arguments of reorder_outputs
% 
%  example:
% 
%     x=rand(3,4); 
% 
%     [cols,rows]=reorder_outputs(@size,[2,1],x); 
% 
%