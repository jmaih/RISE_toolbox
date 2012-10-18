function obj=save_filters(obj)

Fields={'filtered_variables','Expected_filtered_variables',...
    'updated_variables','Expected_updated_variables',...
    'smoothed_variables','Expected_smoothed_variables',...
    'smoothed_shocks','Expected_smoothed_shocks',...
    'smoothed_measurement_errors','Expected_smoothed_measurement_errors',...
    'filtered_probabilities','updated_probabilities',...
    'smoothed_probabilities'};

nobs=numel(obj.varobs(1).value);
reg_names=cell(1,obj.NumberOfRegimes);
for ii=1:obj.NumberOfRegimes
    reg_names{ii}=['regime_',int2str(ii)];
end

for ifield=1:numel(Fields)
    ff=Fields{ifield};
    if isfield(obj.Filters,ff)
        obj.Filters.(ff)=push_filter_to_time_series(obj.Filters.(ff),ff);
    end
end


    function out=push_filter_to_time_series(filtdata,filtname)
        if isempty(filtdata)
            out.filtname=[];
            return
        end
        variables=~isempty(strfind(filtname,'variables'));
        shocks=~isempty(strfind(filtname,'shocks'));
        measerrs=~isempty(strfind(filtname,'measurement_errors'));
        filtered=~isempty(strfind(filtname,'filtered'));
        if filtered
            dateInfo=obj.dates_filtering;
        else
            dateInfo=obj.dates_smoothing;
        end
        [r,c,p]=size(filtdata);
        others='';
        if ismember(p,[nobs,nobs+1]) % this is not very robust...
            filtdata=permute(filtdata,[3,1,2]);
            others=reg_names;
        elseif ismember(c,[nobs,nobs+1])
            filtdata=permute(filtdata,[2,1,3]);
        end
        
        out=struct();
        if variables
            nvars=obj.NumberOfEndogenous(2);
            vnames={obj.varendo.name};
        elseif shocks
            vnames={obj.varexo.name};
            nvars=obj.NumberOfExogenous;
            % add the observations on the exogenous observed variables in
            % order to do things in one go. otherwise there will be a crash
            % here...
            net_nexo=nvars-obj.NumberOfObservables(2);
            nloops=size(filtdata,2)/net_nexo;
            current_reg_nbr=size(filtdata,3);
            tmp=nan(nobs,nvars*nloops,current_reg_nbr);
            % now push in the unobserved shocks and the observed ones (only
            % for the first period)
            % locate the observed shocks
            iobs=ismember(vnames,{obj.varobs_exo.name});
            xxx=vertcat(obj.varobs_exo.value);
            if ~isempty(xxx)
                tmp(:,iobs,:)=repmat(transpose(xxx),[1,1,current_reg_nbr]);
            end
            not_iobs=setdiff(1:nvars,find(iobs));
            index=not_iobs;
            for iloop=2:nloops
                index=[index,not_iobs+nvars];
            end
            tmp(:,index,:)=filtdata;
            filtdata=tmp;
        elseif measerrs
            vnames={obj.varobs.name};
            nvars=obj.NumberOfObservables(1);
        else % probabilities
            vnames=obj.state_names;
            nvars=numel(vnames);
        end
        for ivar=1:nvars
            datta=squeeze(filtdata(:,ivar:nvars:end,:));
            if size(datta,2)>nvars
                others='';
            end
            if isempty(others)
                out.(vnames{ivar})=rise_time_series(dateInfo,datta);
            else
                out.(vnames{ivar})=rise_time_series(dateInfo,datta,others);
            end
        end
    end
end
