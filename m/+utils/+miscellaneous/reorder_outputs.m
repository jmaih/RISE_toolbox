function varargout=reorder_outputs(f,order,varargin)
% reorder_outputs : reorders the outputs of a function
%
% ::
%
%   varargout=reorder_outputs(f,order,varargin)
%
% Args:
%
%    f (function handle): function for which one wants to reorder the
%      outputs
%
%    order (vector): order in which the outputs are to be reordered. Note
%      that all elements in the reordering should be in 1:numel(order)
%
%    varargin (varargin): input arguments to the f function
%
% Returns:
%    :
%
%    - **varargout** [varargout]: output arguments of reorder_outputs
%
% example:
%
%    x=rand(3,4); 
%
%    [cols,rows]=reorder_outputs(@size,[2,1],x); 
%
n=numel(order);

if ~all(ismember(order,1:n))
    
    error('order should contain all elements from 1 to numel(order)')
    
end

[vout{1:n}]=f(varargin{:});

vout=vout(order);

varargout=vout(1:nargout);

end