function output=concatenate_series_from_different_models(dbcell)
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

% dbcell is a cell array of structures of time series objects
% all elements should have the same dimensions else there will be a crash
% returns a structure whose fields are the concatenated time series objects

varList=fieldnames(dbcell{1});
mydates_start=dbcell{1}.(varList{1}).start;
output=struct();
nobj=numel(dbcell);
for ivar=1:numel(varList)
    for imod=1:numel(dbcell)
        datta=double(dbcell{imod}.(varList{ivar}));
        if imod==1 && ivar==1
            data_size=size(datta);
            tank=zeros(data_size(1),nobj,data_size(2));
            mod_names=strcat('model_',cellfun(@num2str,num2cell(1:nobj),'uniformOutput',false));
            reg_names=strcat('regime_',cellfun(@num2str,num2cell(1:data_size(2)),'uniformOutput',false));
        end
        tank(:,imod,:)=datta;
    end
    if data_size(2)>1
        for ireg=1:data_size(2)
            output.(reg_names{ireg}).(varList{ivar})=...
                ts(mydates_start,tank(:,:,ireg),mod_names);
        end
    else
        output.(varList{ivar})=ts(mydates_start,tank,mod_names);
    end
end
end
