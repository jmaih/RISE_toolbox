function ppdata_=confidence_regions(obj,simulation_folder,regions)
% confidence_regions -- computes confidence regions for the posterior
% distribution of parameters
%
% Syntax
% -------
% ::
%
%   ppdata_=confidence_regions(obj)
%
%   ppdata_=confidence_regions(obj,simulation_folder)
%
%   ppdata_=confidence_regions(obj,simulation_folder,regions)
%
% Inputs
% -------
%
% - **obj** [rise|dsge|rfvar|svar]: model object
%
% - **simulation_folder** [[]|char|struct]: location of the simulations. If
% empty, it is assumed that the simulations are saved to disc and are
% located in the address found in obj.folders_paths.simulations. If it is a
% "char", this corresponds to the location of the simulation. Otherwise, if
% it is a struct, then the fields of the struct contain the simulations
%
% - **regions** [[]|scalar|vector]: regions for which the one wants to
% compute the quantiles of the distribution. all elements should be
% strictly greater than 0 and smaller than 100.
%
% Outputs
% --------
%
% - **ppdata_** [struct]: structure whose fields are the names of the
% estimated parameters. After which there are the confidence regions
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

if isempty(obj)
    ppdata_=struct();
    return
end

if nargin<3
    regions=[];
    if nargin<2
        simulation_folder=[];
    end
end

if isempty(regions)
    regions=[80,85,90,95,99];
else
    if ~isa(regions,'double')||...
            ~isvector(regions)||...
            any(regions(:)<=0)||...
            any(regions(:)>=100)
        error('regions must be a scalar or a vector with elements in (0,100)')
    end
end

regions=sort(regions);
regions=regions(:); 
p=100-regions;
p=vec(p.').'; % in this way, the result is also a row vector
nregions=numel(regions);
regions=regions(:)'; 

if isempty(simulation_folder)
    simulation_folder=obj.folders_paths.simulations;
end

is_saved_to_disk=ischar(simulation_folder);
if is_saved_to_disk
    W = what(simulation_folder);
    W=W.mat;
    locs=find(strncmp('chain_',W,6));
    if isempty(locs)
        error([mfilename,':: no simulations found'])
    end
    W=strrep(W(locs),'.mat','');
elseif isstruct(simulation_folder)
    W=fieldnames(simulation_folder);
else
    error('wrong specification of input')
end
number_of_matrices=numel(W);

vnames=cellfun(@(x)parser.param_name_to_valid_param_name(x),...
    {obj.estimation.priors.name},'uniformOutput',false);
tex_names={obj.estimation.priors.tex_name};

% create the data
%----------------
npar=numel(vnames);
ppdata_=struct();
for ipar=1:npar
    all_vals=[];
    for m=1:number_of_matrices
        if is_saved_to_disk
            tmp=load([simulation_folder,filesep,W{m}]);
        else
            tmp=simulation_folder.(W{m});
        end
        Params=tmp.Params(ipar,:);
        all_vals=[all_vals;Params(:)]; %#ok<AGROW>
    end
    ppdata_.(vnames{ipar})=do_one_post(ipar);
end

    function ss=do_one_post(ipar)
        ss=struct();
        ss.tex_name=tex_names{ipar};
        ss.conf_regions=mat2cell(prctile(all_vals,p),1,2*ones(1,nregions));
        ss.percent=regions;
    end

end

