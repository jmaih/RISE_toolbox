function check_factorization(R,Sigma)
% R is the impact matrix for the shocks
% Sigma is the covariance matrix of the errors

if ~isequal(size(R),size(Sigma))
    
    error('wrong format for the impact matrix')
    
end

h=size(R,3);

for ireg=1:h
    
    check_one_factorization(R(:,:,ireg),Sigma(:,:,ireg))

end

    function check_one_factorization(R,Sigma)
        
        RR=R*R';
        
        if any(isnan(RR(:)))||~isreal(RR)
            
            error('factorization returned a bad impact matrix')
            
        end
        
        if max(abs(RR(:)-Sigma(:)))>1e-6
            
            error('impact matrix is not a square root of the covariance of errors')
            
        end
        
    end

end
