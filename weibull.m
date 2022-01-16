function f = weibull(v, k, c)
    
    f = k / c .* (v / c) .^(k-1) .* exp(-(v/c).^k);
    
end