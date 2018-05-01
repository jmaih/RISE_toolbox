function [smpl,fsmpl,accept_rate,start] = mh_sampler(start,nsamples,varargin)
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


start.f0=-start.f0;

[smpl,fsmpl,accept_rate,start] = utils.mcmc.mh_sampling_engine(start,nsamples,...
    'logpdf',@(x)-start.objective(x),...
    'proprnd',@(x)start.drawfun(x,start.c*start.CS),varargin{:});

start.f0=-start.f0;
fsmpl=-fsmpl;