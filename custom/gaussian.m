function ys = gaussian(freqs, mu, hgt, sigma)
%       Gaussian function to use for fitting.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency vector to create gaussian fit for.
%       mu, hgt, sigma : doubles
%           Parameters that define gaussian function (centre frequency,
%           height, and standard deviation).
%
%       Returns
%       -------
%       ys :    1xn array
%       Output values for gaussian function.

    ys = hgt*exp(-(((freqs-mu)./sigma).^2) /2);

end