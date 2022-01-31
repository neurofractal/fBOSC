function err = error_model_constr(params, xVals, yVals, peak_type)
    fitted_vals = 0;
    for set = 1:size(params,1)
        switch (peak_type)
            case 1 % Gaussian
                fitted_vals = fitted_vals + gaussian(xVals, params(set,1), params(set,2), params(set,3));
            case 2 % Cauchy
                fitted_vals = fitted_vals + cauchy(xVals, params(set,1), params(set,2), params(set,3));
        end
    end
    err = sum((yVals - fitted_vals).^2);
end