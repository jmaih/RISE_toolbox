function y = nanvar(x,w,dim)
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

if nargin < 2 || isempty(w), w = 0; end

sz = size(x);
if nargin < 3 || isempty(dim)
    % The output size for [] is a special case when DIM is not given.
    if isequal(x,[]), y = NaN(class(x)); return; end

    % Figure out which dimension sum will work along.
    dim = find(sz ~= 1, 1);
    if isempty(dim), dim = 1; end
elseif dim > length(sz)
    sz(end+1:dim) = 1;
end

% Need to tile the mean of X to center it.
tile = ones(size(sz));
tile(dim) = sz(dim);

if isequal(w,0) || isequal(w,1)
    % Count up non-NaNs.
    n = sum(~isnan(x),dim);

    if w == 0
        % The unbiased estimator: divide by (n-1).  Can't do this when
        % n == 0 or 1, so n==1 => we'll return zeros
        denom = max(n-1, 1);
    else
        % The biased estimator: divide by n.
        denom = n; % n==1 => we'll return zeros
    end
    denom(n==0) = NaN; % Make all NaNs return NaN, without a divideByZero warning

    x0 = x - repmat(utils.stat.nanmean(x, dim), tile);
    y = utils.stat.nansum(abs(x0).^2, dim) ./ denom; % abs guarantees a real result

% Weighted variance
elseif numel(w) ~= sz(dim)
    error('Invalid size of weights');
elseif ~(isvector(w) && all(w(~isnan(w)) >= 0))
    error('Invalid weights');
else
    % Embed W in the right number of dims.  Then replicate it out along the
    % non-working dims to match X's size.
    wresize = ones(size(sz)); wresize(dim) = sz(dim);
    wtile = sz; wtile(dim) = 1;
    w = repmat(reshape(w, wresize), wtile);

    % Count up non-NaNs.
    n = utils.stat.nansum(~isnan(x).*w,dim);

    x0 = x - repmat(utils.stat.nansum(w.*x, dim) ./ n, tile);
    y = utils.stat.nansum(w .* abs(x0).^2, dim) ./ n; % abs guarantees a real result
end