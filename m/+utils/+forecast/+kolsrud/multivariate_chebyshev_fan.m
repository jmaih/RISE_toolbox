function [f,mvcb,c]=multivariate_chebyshev_fan(Y,cr,myvars,start_date)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

c=utils.forecast.kolsrud.chebyshev_distance(Y);

mvcb=utils.forecast.kolsrud.multivariate_chebyshev_box(Y,cr);

[~,T,N,ngam]=size(mvcb);

date_numbers=date2serial(start_date):date2serial(start_date)+T-1;

if ~issorted(cr)
    
    error('cr should be sorted')
    
end

f=struct();

chebyshevbox2fan();

    function chebyshevbox2fan()
        dataL=zeros(T,ngam);
        
        dataR=zeros(T,ngam);
        
        mean_=squeeze(mean(Y,1));
        
        for ivar=1:N
            
            for ig=1:ngam
                
                cut=permute(mvcb(:,:,ivar,ig),[2,1]);
                
                dataL(:,end-ig+1)=cut(:,1);
                
                dataR(:,ig)=cut(:,2);
                
            end
            
            f.(myvars{ivar})=struct('quantiles',[dataL,dataR],...
                'ci',cr,'date_numbers',date_numbers,...
                'mean',mean_(:,ivar));
            
        end
    end

end

