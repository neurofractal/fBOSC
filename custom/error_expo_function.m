function err = error_expo_function(params,xs,ys)
    ym = expo_function(xs,params);
    err = sum((ys - ym).^2);
end