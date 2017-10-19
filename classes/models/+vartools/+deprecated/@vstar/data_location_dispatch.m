function obj=data_location_dispatch(obj)

% prevent crashes...
vlist=strrep(obj.observables.name,'{0}','');

d=struct();

tnames=threshold_names();

d.endo_id=locate_variables(obj.endogenous.name,vlist);

d.det_id=locate_variables(obj.deterministic.name,vlist);

d.thresh_id=locate_variables(tnames,vlist);

obj.variables_locations_in_data=d;

    function tnames=threshold_names()
        
        nthresh=numel(obj.thresholds);
        
        tnames=cell(1,nthresh);
        
        for ithresh=1:nthresh
            
            thresh=obj.thresholds(ithresh);
            
            vname=thresh.name;
            
            if ~strcmp(thresh.lag,'0')
                
                vname=[vname,'{',thresh.lag,'}']; %#ok<AGROW>
                
            end
            
            tnames{ithresh}=vname;
            
        end
        
    end

end