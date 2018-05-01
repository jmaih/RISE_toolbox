function out=empirical_moments(xx,lb,ub,probs,kernel,number_of_grid_points)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:


if nargin<6
	number_of_grid_points=[];
	if nargin<5
		kernel=[];
		if nargin<4
			probs=[];
			if nargin<3
				ub=[];
				if nargin<2
					lb=[];
				end
			end
		end
	end
end

if isempty(kernel)
	kernel='normal';
end
if isempty(lb)
	lb=min(xx,[],2);
end
if isempty(ub)
	ub=max(xx,[],2);
end
if isempty(number_of_grid_points)
    number_of_grid_points = 2^9;      % 2^9 = 512 !... Must be a power of two.
end
if isempty(probs)
    probs = 0.1:.1:.9;
end

[npar,number_of_draws]=size(xx);
if numel(lb)~=npar||numel(ub)~=npar
	error('number of elements in bounds inconsistent with the number of parameters as given by xx')
end

xx=sort(xx,2);
out_mean = mean(xx,2);
out_median = median(xx,2);
out_variance = var(xx,[],2);
out_min = min(xx,[],2);
out_max = max(xx,[],2);
prob_index=round(probs*number_of_draws);
out_quantiles = xx(:,prob_index);

out=struct();
for ipar=1:npar
	out(ipar).mean=out_mean(ipar);
	out(ipar).min=out_min(ipar);
	out(ipar).max=out_max(ipar);
	out(ipar).median=out_median(ipar);
	out(ipar).variance=out_variance(ipar);
	out(ipar).quantiles=out_quantiles(ipar,:);
	out(ipar).prob_index=prob_index;
	out(ipar).probs=probs;
	[out(ipar).fdens,out(ipar).xdens]=distributions.kernel_density(xx(ipar,:),lb(ipar),ub(ipar),kernel,number_of_grid_points);
	[out(ipar).cdf]=distributions.empirical_cdf(xx(ipar,:),lb(ipar),ub(ipar),number_of_grid_points);
end
