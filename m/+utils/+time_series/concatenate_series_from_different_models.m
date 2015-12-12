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
    is_good=true;
    for imod=1:numel(dbcell)
        is_good=is_good && isfield(dbcell{imod},varList{ivar});
        if ~is_good
            continue
        end
        this=dbcell{imod}.(varList{ivar});
        dn=this.date_numbers;
        datta=double(this);
        if imod==1 && ivar==1
            data_size=size(datta);
            tank=nan(data_size(1),nobj,data_size(2));
            mod_names=strcat('model_',cellfun(@num2str,num2cell(1:nobj),'uniformOutput',false));
            reg_names=strcat('regime_',cellfun(@num2str,num2cell(1:data_size(2)),'uniformOutput',false));
            dn0=dn;
            big_start=dn0(1);
            big_end=dn0(end);
            tank(:,imod,:)=datta;
        else
            before=big_start-dn(1);
            if before>0
                big_start=dn(1);
                tank=cat(1,nan(before,nobj,data_size(2)),tank);
            end
            after=dn(end)-big_end;
            if after>0
                big_end=dn(end);
                tank=cat(1,tank,nan(after,nobj,data_size(2)));
            end
            if before||after
                dn0=big_start:big_end;
            end
            % better start search from the beginning
            start=find(dn(1)==dn0,1,'first');
            % better start search from the end
            finish=find(dn(end)==dn0,1,'last');
            tank(start:finish,imod,:)=datta;
        end
    end
    if ~is_good
        disp(['skipping "', varList{ivar},'" as it is not included in all models'])
        continue
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
