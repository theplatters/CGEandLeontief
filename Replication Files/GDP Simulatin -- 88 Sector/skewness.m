function sk = skewness(x)
sk = (sum((x-mean(x)).^3)./length(x)) ./ (var(x,1).^1.5);
end