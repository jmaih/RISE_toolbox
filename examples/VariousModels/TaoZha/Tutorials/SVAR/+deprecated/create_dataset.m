function [db,varlist,fh]=create_dataset(scale,do_plot)
% create_dataset -- creates time series for all SVAR models
%
% ::
%
%
%   [db,varlist,fh]=create_dataset()
%
%   [db,varlist,fh]=create_dataset(scale)
%
%   [db,varlist,fh]=create_dataset(scale,do_plot)
%
% Args:
%
%    - **scale** [empty|numeric|{1}]: factor by which to multiply the data
%
%    - **do_plot** [empty|true|{false}]: if true, the data are plotted.
%
% Returns:
%    :
%
%    - **db** [struct]: structure with fields as the names of the endogenous
%    variables.
%
%    - **varlist** [cellstr]: list of the endogenous variables and their
%    description
%
%    - **fh** [empty|handle]: if **do_plot** is set to true, **fh** is the
%    handle to the plotted figure. Else, it is empty
%
% Note:
%
% Example:
%
%    See also:

if nargin<2
    
    do_plot=false;
    
    if nargin<1
        
        scale=1;
        
    end
    
end

if isempty(do_plot),do_plot=false; end

if isempty(scale),scale=1; end

if scale<=0
    
    error('scale must be strictly positive')
    
end

yrBin  = 1954;   % beginning year

qmBin  = 3;      % begining quarter or month

% names of the variables
%------------------------
varlist=struct();

varlist.FFR='Feds Funds Rate';

varlist.pi='Inflation';

varlist.ygap='Output gap';

% RISE uses "" for the description of the model variables
%--------------------------------------------------------

rawdb=load('dataraw_allvars.mat');

start_date=sprintf('%0.0dQ%0.0d',yrBin,qmBin);

% create the data as a page
%---------------------------
db=ts(start_date,scale*rawdb.xdd,fieldnames(varlist));

% separate the various variables
%-------------------------------
db=pages2struct(db);

fh=[];

if do_plot
    
    fh=figure('name','Variables in the VAR');
    
    fields=fieldnames(varlist);
    
    for iplot=1:3
        
        vname=fields{iplot};
        
        subplot(3,1,iplot)
        
        plot(db.(vname),'linewidth',2)
        
        mytitle=strrep([varlist.(vname),'(',vname,')'],'"','');
        
        title(mytitle)
        
    end
    
end

end