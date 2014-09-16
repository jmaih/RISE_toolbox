function [h,legend_]=prior_posterior(ss,varargin)
% initialize the legend items
%----------------------------
legend_={};

% plot posterior
%---------------
plot(ss.x_kdens,ss.f_kdens,'LineStyle','-','Color','b',varargin{:})
legend_=[legend_,'post density'];
hold on

% plot prior
%-----------
plot(ss.x_prior,ss.f_prior,'LineStyle','-','Color','green',varargin{:})
legend_=[legend_,'prior density'];

% vertical line at posterior mean
%--------------------------------
[x_mean,position]=utils.miscellaneous.find_nearest(ss.x_kdens,ss.mean_sim);
plot([x_mean x_mean], [0 ss.f_kdens(position)],...
    'LineStyle',':','Color','black',varargin{:} )
legend_=[legend_,'mean'];

% plot vertical line at maximization mode
%---------------------------------------
[x_mode,position]=utils.miscellaneous.find_nearest(ss.x_kdens,ss.post_mode);
plot([x_mode x_mode],[0 ss.f_kdens(position)],...
    'LineStyle',':','Color','green',varargin{:} )
legend_=[legend_,'max-mode'];

% plot vertical line at simulation mode
%---------------------------------------
[x_post_mode_sim,position]=utils.miscellaneous.find_nearest(ss.x_kdens,ss.post_mode_sim);
plot([x_post_mode_sim x_post_mode_sim],[0 ss.f_kdens(position)],...
    'LineStyle',':','Color','red',varargin{:} )
legend_=[legend_,'sim-mode'];

hold off

% size of the plot
%-----------------
xlow=min(ss.min_sim,min(ss.x_prior));
xhigh=max(ss.max_sim,max(ss.x_prior));
axis([xlow xhigh 0 1.4*max(ss.f_kdens)]);

h=gca();
end
