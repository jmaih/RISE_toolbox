function [db,varlist,fh]=create_dataset(do_plot)
% create_dataset -- creates time series for all SVAR models
%
% Syntax
% -------
% ::
%
%   [db,varlist,fh]=create_dataset()
%
%   [db,varlist,fh]=create_dataset(do_plot)
%
% Inputs
% -------
%
% - **do_plot** [empty|true|{false}]: if true, the data are plotted.
%
% Outputs
% --------
%
% - **db** [struct]: structure with fields as the names of the endogenous
% variables.
%
% - **varlist** [cellstr]: list of the endogenous variables
%
% - **fh** [empty|handle]: if **do_plot** is set to true, **fh** is the
% handle to the plotted figure. Else, it is empty
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

if nargin==0
    do_plot=false;
end

yrBin  = 1954;   % beginning year

qmBin  = 3;      % begining quarter or month

% names of the variables
%------------------------
varlist   = {'FFR','pi','ygap'};

rawdb=load('dataraw_allvars.mat');

start_date=sprintf('%0.0dQ%0.0d',yrBin,qmBin);

% create the data as a page
%---------------------------
db=ts(start_date,rawdb.xdd,varlist);

% separate the various variables
%-------------------------------
db=pages2struct(db);

fh=[];
if do_plot
    fh=figure('name','Variables in the VAR');
    for iplot=1:3
        vname=varlist{iplot};
        subplot(3,1,iplot)
        plot(db.(vname),'linewidth',2)
        title(vname)
    end
end

end