function err = error_expo_nk_function(params,xs,ys)
    ym = -log10(xs.^params(2)) + params(1);
    err = sum((ys - ym).^2);
end
