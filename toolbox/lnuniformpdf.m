function logpdf = lnuniformpdf(x, prior)
%
%


logpdf = log(abs(prior.constraint(2) - prior.constraint(1)));   