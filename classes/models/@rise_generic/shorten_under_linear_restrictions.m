function varargout=shorten_under_linear_restrictions(obj,varargin)
% shorten_under_linear_restrictions -- shorten the size of vectors under
% linear restrictions
%
% Syntax
% -------
% ::
%
%   varargout=shorten_under_linear_restrictions(obj,varargin)
%
% Inputs
% -------
%
% - **obj** [rise|dsge|rfvar|svar]: model object
%
% - **varargin** [vector|2-column matrix|cell]: elements to shorten. In
% case of a two-column matrix, it is assumed that the first column is the
% lower bound of the parameters and the second column is the upper bound.
% Some checks have to be made after transformation in order to insure that
% no element in the transformed lower bound exceeds its upper bound
% counterpart. In case of a cell, it is assumed that we are dealing with a
% covariance matrix.
%
% Outputs
% --------
%
% - **varargout** [vector|2-column matrix]: transformed inputs
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    varargout={struct()};
    return
end

% shorten everything in the presence of linear restrictions
%-----------------------------------------------------------
n=length(varargin);
varargout=cell(1,n);
linear_restricts=obj(1).linear_restrictions_data;
for ii=1:n
    if iscell(varargin{ii})
        % covariance matrix
        vcov=varargin{ii}{1};
        if isempty(vcov)
            continue
        end
        if size(vcov,1)~=size(vcov,2)
            error('expected a covariance matrix... sizes do not match')
        end
        % the covariance matrix has a special call to the function
         varargout{ii}=linear_restricts.a2tilde_func(vcov,true);
    else
        [~,ncols]=size(varargin{ii});
        if ncols==0
            continue
        elseif ncols==1
            varargout{ii}=linear_restricts.a2tilde_func(varargin{ii});
        elseif ncols==2 % bounds
            lb=linear_restricts.a2tilde_func(varargin{ii}(:,1));
            ub=linear_restricts.a2tilde_func(varargin{ii}(:,2));
            bad=ub<lb;
            if any(bad)
                tmp=lb;
                lb(bad)=ub(bad);
                ub(bad)=tmp(bad);
            end
            varargout{ii}=[lb,ub];
        else
            error('unrecognized input')
        end
    end
end

end