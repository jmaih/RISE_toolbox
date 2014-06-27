function out=initialize_solution_or_structure(type,nreg)

out=struct();
if nargin==0
    return
end
switch type
    case {'sol','solution'}
        prototype=cell(1,nreg);
        mylist={'definitions','ss','bgp',...
            'm_x',... % formerly T
            'm_e',... % formerly R
            'm_sig',...
            'm_x_x','m_x_e','m_x_sig','m_e_e','m_e_sig','m_sig_sig','theta_hat','Q','H'};
    case {'system','structure'}
        prototype=cell(nreg);
        mylist={'Gp','Gc','Gm','Ge','Gt','Gpp','Gpc','Gpm','Gpe','Gpt','Gcp',...
            'Gcc','Gcm','Gce','Gct','Gmp','Gmc','Gmm','Gme','Gmt','Gep',...
            'Gec','Gem','Gee','Get','Gtp','Gtc','Gtm','Gte','Gtt'};
    otherwise
        error(['unknown type: ',type])
end

for ilist=1:numel(mylist)
    d=mylist{ilist};
    if strcmp(d,'Q')
        out.(d)=sparse(nreg,nreg);
    else
        out.(d)=prototype;
    end
end

if ismember(type,{'system','structure'})
    out.planner=struct('discount',[],'commitment',[],'objective',[],'weights',[]);
end

