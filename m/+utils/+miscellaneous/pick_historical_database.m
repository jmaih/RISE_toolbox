function histdb=pick_historical_database(m,type,hist_end_date)
% pick_historical_database pick a historical database in the filters
%
% Syntax
% -------
% ::
%
%   histdb=pick_historical_database(m)
%
%   histdb=pick_historical_database(m,type)
%
%   histdb=pick_historical_database(m,type,hist_end_date)
%
% Inputs
% -------
%
% - **m** [rise|dsge] : model object
%
% - **type** [{'mean'}|'random'] : type of database. If 'mean', the
% database with expected values is returned. If 'random', a database is
% drawn randomly using the probabilitiy distribution at the end of history.
%
% - **hist_end_date** [char|{sample end date}] : desired date for the end
% of history
%
% Outputs
% --------
%
% - **histdb** [struct] : structure with the time series of the endogenous
% variables, the regime probabilities
%
% More About
% ------------
%
% - if the number of regimes is 1, "type" and "hist_end_date" are
% superfluous.
%
% - if a database is picked at random, the regime probabilities are changed
% so as to inform the forecasting procedure later on about what regime the
% data belong to. For instance, suppose we have 3 regimes and that the
% second regime is picked based on the distribution [.5,.2,.3] at the
% desired (or default) end date of history. Once the regime is picked,
% those probabilities are modified to [0,1,0].
%
% Examples
% ---------
%
% See also: 

if nargin<3
    hist_end_date=[];
    if nargin<2
        type=[];
    end
end

if ~isa(m,'dsge')
    error('First input must be a dsge or a rise object')
end

if isempty(hist_end_date)
    hist_end_date=m.options.estim_end_date;
end
regimes_number=m.markov_chains.regimes_number;

if regimes_number==1||isempty(type)
    type='mean';
end

if ~any(strcmpi(type,{'mean','random'}))
    error('second input must be "random" or "mean"')
end

switch lower(type)
    case 'mean'
        histdb=m.filtering.Expected_updated_variables;
        histdb=utils.miscellaneous.mergestructures(histdb,...
            m.filtering.updated_regime_probabilities);
    case 'random'
        reg_probs0=struct2pages(m.filtering.updated_regime_probabilities);
        reg_probs=cumsum(double(reg_probs0(hist_end_date)));
        reg_probs=[0,reg_probs(:).'];
        pick=find(reg_probs>rand,1,'first')-1;
        histnames=m.endogenous.name;
        histdb=struct();
        reg=sprintf('regime_%0.0f',pick);
        for iname=1:numel(histnames)
            name=histnames{iname};
            histdb.(name)=m.filtering.updated_variables.(name)(reg);
        end
        hist_probs=zeros(1,regimes_number);
        hist_probs(pick)=1;
        reg_probs0(hist_end_date,:)=hist_probs;
        reg_probs0=pages2struct(reg_probs0);
        histdb=utils.miscellaneous.mergestructures(histdb,reg_probs0);
end


end