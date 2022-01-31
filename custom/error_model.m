function err = error_model(params, xVals, yVals, peak_type, guess, guess_weight)
    fitted_vals = 0;
    weak = 1E2;
    strong = 1E7;
    for set = 1:size(params,1)
        switch (peak_type)
            case 1 % Gaussian
                fitted_vals = fitted_vals + gaussian(xVals, params(set,1), params(set,2), params(set,3));
            case 2 % Cauchy
                fitted_vals = fitted_vals + cauchy(xVals, params(set,1), params(set,2), params(set,3));
        end
    end
    switch guess_weight
        case 'none'
            err = sum((yVals - fitted_vals).^2);
        case 'weak' % Add small weight to deviations from guess m and amp
            err = sum((yVals - fitted_vals).^2) + ...
                 weak*sum((params(:,1)-guess(:,1)).^2) + ...
                 weak*sum((params(:,2)-guess(:,2)).^2);
        case 'strong' % Add large weight to deviations from guess m and amp
            err = sum((yVals - fitted_vals).^2) + ...
                 strong*sum((params(:,1)-guess(:,1)).^2) + ...
                 strong*sum((params(:,2)-guess(:,2)).^2);
    end
end