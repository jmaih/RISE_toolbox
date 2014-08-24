function x1=build_parameter_vector(vdata,a_post,vcov)

persistent same npar p1 p2 nan_locs

initialize=isempty(npar);

if initialize
    npar=numel(vdata.orig_estim_names);
end
x1=nan(npar,1);
x1(vdata.estim_locs)=a_post;
if initialize
    nan_locs=find(isnan(x1));
    pp=regexp(vdata.orig_estim_names(nan_locs),...
        '\w+(?<first>\d+)_(?<second>\d+)','names');
    pp=[pp{:}];
    p1=cell2mat(cellfun(@(x)str2double(x),{pp.first},'uniformOutput',false));
    p2=cell2mat(cellfun(@(x)str2double(x),{pp.second},'uniformOutput',false));
    same=p1==p2;
end

[SIG,OMG]=vartools.covariance_decomposition(vcov);

pvals=nan(size(same));
pvals(same)=SIG(p1(same));
pvals(~same)=diag(OMG(p1(~same),p2(~same)));

x1(nan_locs)=pvals;

end

% %     for ipar=1:numel(pp)
% %         first=str2double(pp(ipar).first);
% %         second=str2double(pp(ipar).second);
% %         if first==second
% %             x1(nan_locs(ipar))=sdev(first);
% %         else
% %             x1(nan_locs(ipar))=vcov(first,second)/(sdev(first)*sdev(second));
% %         end
% %     end
% %     for irow=1:nvars
% %         sig_name=sprintf('sig_%0.0f_%0.0f',irow,irow);
% %         sig_pos=strcmp(sig_name,vdata.orig_estim_names);
% %         x1(sig_pos)=sdev(irow);
% %         for icol=1:irow-1
% %             corr_name=sprintf('omg_%0.0f_%0.0f',irow,icol);
% %             corr_pos=strcmp(corr_name,vdata.orig_estim_names);
% %             x1(corr_pos)=vcov(irow,icol)/(sdev(irow)*sdev(icol));
% %         end
% %     end
% %     keyboard
