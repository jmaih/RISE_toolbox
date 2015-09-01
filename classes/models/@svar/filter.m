function [obj,LogLik,Incr,retcode]=filter(obj,varargin)
% FILTER -- computes historical probabilities and smoothed shocks for VARs
%
% Syntax
% -------
% ::
%
%   [obj,LogLik,Incr,retcode]=filter(obj,varargin)
%
% Inputs
% -------
%
% - **obj** [svar|rfvar]: model object
%
% - **varargin** []: pairwise options
%
% Outputs
% --------
%
% - **obj** [svar|rfvar]: model object with the filters
%
% - **LogLik** [numeric]: log likelihood
%
% - **Incr** [vector]: contributions to the likelihood for each time t
%
% - **retcode** [integer]: 0 if no problem encountered
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: dsge/filter

if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    obj=struct();
    return
end

nobj=numel(obj);
Incr=[];
if nobj>1
    retcode=nan(1,nobj);
    LogLik=nan(1,nobj);
    for iobj=1:nobj
        [obj(iobj),LogLik(iobj),~,retcode(iobj)]=filter(obj(iobj),varargin{:});
    end
    return
end

if ~isempty(varargin)
    obj=set(obj,varargin{:});
end

if ~obj.data_are_loaded
    obj=obj.load_data;
end

[LogLik,Incr,retcode,obj]=vartools.var_likelihood([],obj);

end