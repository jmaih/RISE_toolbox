function check_errors(mydata,nlags)

if size(mydata,1)<0
    
    error('need at least 1 variable in your data')
    
end

if nlags<1
    
    error('need at least 1 lag; set nlags higher')
    
end

if nlags>size(mydata,2)
    
    error('need more data points than lags')
    
end

if any(isnan(mydata(:)))
    
    error('no missing observations allowed')
    
end

end
