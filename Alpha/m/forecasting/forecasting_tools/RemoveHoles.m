function [rstar,varargout]=RemoveHoles(rstar,varargin)
holes=isnan(rstar);
rstar(holes,:)=[];
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

