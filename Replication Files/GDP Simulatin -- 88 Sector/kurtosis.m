function kurtss = kurtosis(x)
kurtss =(sum((x-mean(x)).^4)./length(x)) ./ (var(x,1).^2);
end