% clear all
% load('Yn.mat')
y = Yn(:,2); 
% % Likelihood function for Gaussian distribution
% likelihood_function = @(params) -sum(log(normpdf(y, params(1), params(2))));
% initial_params = [0, 1000];
% 
% % Use fminunc to estimate the parameters that maximize the likelihood
% options = optimset('fminunc');
% options.Display = 'off';
% estimated_params = fminunc(@(params) likelihood_function(params), initial_params, options);

% Estimated parameters
% mu_estimated = estimated_params(1);
% sigma_estimated = estimated_params(2);
% 
% % Calculate E(Y^2)
% e_y_squared = mu_estimated^2 + sigma_estimated^2;

% Calculate the Hessian matrix (second partial derivatives)
% syms sigma sigma_epsilon t_i r_i;  % Define symbolic parameters
% likelihood = 1 / (sigma_epsilon * sqrt(2 * pi)) * exp(-(r_i^2) / (2 * sigma_epsilon^2));
% log_likelihood = sum(log(subs(likelihood, {r_i, sigma_epsilon}, {y,sigma})));
% d2log_likelihood_dsigma = diff(log_likelihood, sigma, 2);
% % d2log_likelihood_dmu = diff(log_likelihood, mu, 2);
% d2log_likelihood_dsigma_dmu = diff(diff(log_likelihood, sigma), mu);
% Fisher_Information_Matrix = [d2log_likelihood_dmu,d2log_likelihood_dsigma_dmu; ...
%                             d2log_likelihood_dsigma_dmu, d2log_likelihood_dsigma ];
% hessian_matrix=subs(Fisher_Information_Matrix, [mu,sigma], [mu_estimated, sigma_estimated]);
% 
% variance_e_y_squared = e_y_squared * hessian_matrix(1, 1) + ...
%                       2 * e_y_squared * hessian_matrix(1, 2) + ...
%                       hessian_matrix(2, 2)^2;
% variance_e_y_squared=-1/(d2log_likelihood_dsigma)
% variance_e_y_squared=double(subs(variance_e_y_squared,[sigma],[sigma_estimated]))
% disp(['Estimated E(Y^2): ' num2str(e_y_squared)]);
% disp(['Variance of E(Y^2): ' num2str(variance_e_y_squared)]);

Mdl = garch('GARCHLags',1,'ARCHLags',1,'Offset',NaN);
[EstMdl,EstParamCov] = estimate(Mdl, y);
disp('Estimated GARCH(1,1) Parameters:');
disp(EstMdl);
% Extract parameter estimates for alpha and beta
estimated_nui_Alpha = EstMdl.ARCH{1};
estimated_nui_Beta = EstMdl.GARCH{1};
estimated_nui_parameters = [estimated_nui_Alpha, estimated_nui_Beta];
%estimated_nui_omega=e_y_squared*(1-estimated_nui_Alpha-estimated_nui_Beta )
estimated_nui_omega=EstMdl.Constant
vv=EstParamCov(1,1)
