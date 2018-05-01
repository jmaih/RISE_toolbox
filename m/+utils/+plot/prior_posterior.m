function [h,legend_,tex_name]=prior_posterior(ss,varargin)
% prior_posterior -- plots prior and, if information is available,
% posterior densities
%
% ::
%
%
%   [h,legend_]=prior_posterior(ss,varargin)
%
% Args:
%
%    - **ss** [struct]: structure containing the relevant information to plot
%    and more specifically
%      - **x_kdens** [vector]: x-axis values for the posterior density
%      - **f_kdens** [vector]: y-axis values for the posterior density
%      - **x_prior** [vector]: x-axis values for the prior density
%      - **f_prior** [vector]: y-axis values for the prior density
%      - **mean_sim** [scalar]: mean of the posterior simulation
%      - **post_mode** [scalar]: value at the posterior maximization mode
%      - **post_mode_sim** [scalar]: value at the posterior simulation mode
%      - **tex_name** [char]: name of the parameter
%
%    - **varargin** [pairwise arguments]: standard plotting arguments for
%    matlab
%
% Returns:
%    :
%
%    - **h** [handle]: handle for the plot
%
%    - **legend_** [cellstr]: names of the lines in the plot
%
%    - **tex_name** [char]: name of the parameter
%
% Note:
%
%    - Only the prior density is plotted if no posterior information is
%    available.
%
% Example:
%
%    See also:

% initialize the legend items
%----------------------------
legend_={};

top=0;
bottom=inf;
is_prior=isfield(ss,'x_prior');
is_posterior=isfield(ss,'x_kdens');
if is_posterior
    % plot posterior
    %---------------
    plot(ss.x_kdens,ss.f_kdens,'LineStyle','-','Color','b',varargin{:})
    legend_=[legend_,'post density'];
    hold on
    top=max(top,max(ss.f_kdens));
    bottom=min(bottom,min(ss.f_kdens));
end

% plot prior
%-----------
if is_prior
    plot(ss.x_prior,ss.f_prior,'LineStyle','-','Color','green',varargin{:})
    legend_=[legend_,'prior density'];
    top=max(top,max(ss.f_prior));
    bottom=min(bottom,min(ss.f_prior));
end

if is_posterior
    % vertical line at posterior mean
    %--------------------------------
    [x_mean,position]=utils.miscellaneous.find_nearest(ss.x_kdens,ss.mean_sim);
    plot([x_mean x_mean], [0 ss.f_kdens(position)],...
        'LineStyle',':','Color','black',varargin{:} )
    legend_=[legend_,'mean'];
    
    % plot vertical line at maximization mode
    %---------------------------------------
    if isfield(ss,'post_mode')
        [x_mode,position]=utils.miscellaneous.find_nearest(ss.x_kdens,ss.post_mode);
        plot([x_mode x_mode],[0 ss.f_kdens(position)],...
            'LineStyle',':','Color','green',varargin{:} )
        legend_=[legend_,'max-mode'];
    end
    
    % plot vertical line at simulation mode
    %---------------------------------------
    [x_post_mode_sim,position]=utils.miscellaneous.find_nearest(ss.x_kdens,ss.post_mode_sim);
    plot([x_post_mode_sim x_post_mode_sim],[0 ss.f_kdens(position)],...
        'LineStyle',':','Color','red',varargin{:} )
    legend_=[legend_,'sim-mode'];
    
    hold off
end

% this title may be over-written by the caller!!!
%------------------------------------------------
title(ss.tex_name)

tex_name=ss.tex_name;

% size of the plot
%-----------------

if abs(top-bottom)<sqrt(eps)
    
    top=inf;
    
end

axis([ss.x_min ss.x_max bottom 1*top]);

h=gca();
end
