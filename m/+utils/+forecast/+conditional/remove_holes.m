function [rstar,varargout]=remove_holes(rstar,varargin)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

holes=isnan(rstar);
rstar(holes,:)=[];
varargout=varargin;
for j=1:length(varargin)
    V=varargin{j};
    % check whether V is a covariance matrix
    if iscell(V) % then it is a covariance matrix
        V=[V{:}];
        if ~isempty(V)
            varargout{j}=V(~holes,~holes);
        else
            varargout{j}=V;
        end
    else
        V(holes,:)=[];
        varargout{j}=V;
    end
end

