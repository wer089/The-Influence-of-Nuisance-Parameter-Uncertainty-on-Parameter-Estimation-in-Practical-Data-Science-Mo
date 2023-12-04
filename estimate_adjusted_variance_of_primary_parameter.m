function estimate_adjusted_variance_of_primary_parameter(returns,setOmega,estimated_nui_omega,variance_e_y_squared,estimated_nui_Alpha,estimated_nui_Beta)
    % Define the log_likelihood function for nuisane parameter
    % the reference of log_likehood function comes from https://www.jstor.org/stable/3318818
    syms omega a b sigma_epsilon t_i r_i;
    likelihood_omega = r_i^2/sigma_epsilon^2+log(sigma_epsilon^2);
    log_likelihood_omega = mean(subs(likelihood_omega, {r_i, sigma_epsilon}, {returns, omega/(1 - a - b)}));
    d2log_likelihood_domega = diff(log_likelihood_omega, omega, 1);
    
    %% Define the log_likelihood function for primary parameter
    %reference of sigma^2:https://math.berkeley.edu/~btw/thesis4.pdf
    syms omega a b sigma_epsilon t_i r_i; 
    likelihood = (1 / sqrt(2 * pi * sigma_epsilon)) * exp(-(r_i^2) / (2 * sigma_epsilon));
    log_likelihood = mean(log(subs(likelihood, {r_i, sigma_epsilon}, {returns, omega/(1 - a - b)})));
    d2log_likelihood_dalpha_dalpha = diff(log_likelihood, a, 2);
    d2log_likelihood_dbeta_dbeta = diff(log_likelihood, b, 2);
    d2log_likelihood_dalpha_dbeta = diff(diff(log_likelihood, a),b);
    d2log_likelihood_dalpha_domega = diff(diff(log_likelihood, a),omega);
    d2log_likelihood_dbeta_domega = diff(diff(log_likelihood, b),omega);
    Hessian_Matrix = [d2log_likelihood_dalpha_dalpha, 0; ...
                               0, d2log_likelihood_dbeta_dbeta];
    %% 0.Estimate the GARCH(1,1) model with a fixed constant (omega)
    estimatedOmega = estimated_nui_omega;
    %estimatedOmega=0.02
    Mdl = garch('Constant',estimatedOmega,'GARCHLags',1,'ARCHLags',1,'Offset',NaN);
    [EstMdl,EstParamCov] = estimate(Mdl, returns);
    disp('Estimated GARCH(1,1) Parameters:');
    disp(EstMdl);
    % Extract parameter estimates for alpha and beta
    estimatedAlpha = EstMdl.ARCH{1};
    estimatedBeta = EstMdl.GARCH{1};
    estimated_parameters = [estimatedAlpha, estimatedBeta];;
    cov0 = EstParamCov(2:3,2:3);
    
    %% 1.Estimate the GARCH(1,1) model with a fixed constant (omega)
    if setOmega==NaN
        Mdl = garch('GARCHLags',1,'ARCHLags',1,'Offset',NaN);
    else
        Mdl = garch('Constant',setOmega,'GARCHLags',1,'ARCHLags',1,'Offset',NaN);
    [EstMdl,EstParamCov] = estimate(Mdl, returns);
    disp('Estimated GARCH(1,1) Parameters:');
    disp(EstMdl);
    % Extract parameter estimates for alpha and beta
    estimatedAlpha = EstMdl.ARCH{1};
    estimatedBeta = EstMdl.GARCH{1};
    estimated_parameters = [estimatedAlpha, estimatedBeta];
    cov = EstParamCov(2:3,2:3);
    
    %% adjusted cov for 1
    %estimatedOmega=sigma_estimated*(1-estimatedAlpha-estimatedBeta);
    c=double(subs([d2log_likelihood_dalpha_domega;d2log_likelihood_dbeta_domega],[a, b,omega], [estimatedAlpha, estimatedBeta, estimatedOmega]));
    D_sym=cov0*[d2log_likelihood_dalpha_domega;d2log_likelihood_dbeta_domega];
    %D_sym=cov0*[d2log_likelihood_domega;d2log_likelihood_domega];
    D=double(subs(D_sym,[a, b,omega], [estimatedAlpha, estimatedBeta,estimatedOmega]));
    v=variance_e_y_squared*(1-estimated_nui_Alpha-estimated_nui_Beta)^2;
    adjust=D*v*transpose(D);
    new_cov=cov0+D*v*transpose(D);

    %% Plot the covariance ellipse for 1
    critical_value = chi2inv(0.90, 2);
    figure;
    hold on
    plot_ellipse(critical_value, cov, estimated_parameters,'r');
    plot_ellipse(critical_value, cov0, estimated_parameters,'black');
    plot_ellipse(critical_value, new_cov, estimated_parameters,'b');
    hold on
    scatter(0.25, 0.7, 'green', 'filled' )
    axis equal;
    xlabel('Alpha');
    ylabel('Beta');
    legend('unadjusted', 'true value', 'adjusted');
    title({['90% Confidence Sigma Ellipse for Alpha and Beta:'];['omega=',num2str(setOmega),', true omega=',num2str(estimatedOmega)]});
    grid on;
    hold off
end
    