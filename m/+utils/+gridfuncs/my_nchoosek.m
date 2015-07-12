function c = my_nchoosek(n,k)

    if k > n
        error('k cannot be bigger than n');
    end
    
    if k > n/2   % use smaller k if available
        k = n-k;
    end
    
    if k <= 1
        c = n^k;
    else
        nums = (n-k+1):n;
        dens = 1:k;
        nums = nums./dens;
        c = round(prod(nums));
    end
    
end
