function [eepd,ep,fh]=setup_endogenous_priors_engine(prior_trunc,fh)
% INTERNAL FUNCTION
%

eepd=[];

ep=[];

if isempty(fh)

    return

end

if ~isa(fh,'function_handle')

    error('estim_endogenous_priors must be a function handle')

end

test=fh();

test_start=false; % do not test whether the start value is within bounds

ep=[];

for ii=1:numel(test.priors)

    % add an initial value: it doesn't matter but it is necessary in order
    % to avoid a crash
    test1=[{0.5},test.priors{ii}];

    blk=utils.prior.cell2block(test1,['endo_prior_',int2str(ii)],'const',1,[]);

    ep=utils.prior.prior_setting_engine(ep,blk,1,prior_trunc,test_start);

end

eepd=utils.prior.load_priors(struct(),ep,[]);

end