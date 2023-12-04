% clear all
% load('Yn.mat')
returns = Yn(:,2);
%vv=3000000000
% Define the log_likelihood function for nuisane parameter
% the reference of log_likehood function comes from https://www.jstor.org/stable/3318818
syms omega alpha beta sigma_epsilon t_i r_i;
likelihood_omega = r_i^2/sigma_epsilon+log(sigma_epsilon);
log_likelihood_omega = mean(subs(likelihood_omega, {r_i, sigma_epsilon}, {returns, omega/(1 - alpha - beta)}));
d2log_likelihood_domega = diff(log_likelihood_omega, omega, 1);

%% Define the log_likelihood function for primary parameter
%reference of sigma^2:https://math.berkeley.edu/~btw/thesis4.pdf
syms omega alpha beta sigma_epsilon t_i r_i; 
likelihood = (1 / sqrt(2 * pi * sigma_epsilon)) * exp(-(r_i^2) / (2 * sigma_epsilon));
log_likelihood = mean(log(subs(likelihood, {r_i, sigma_epsilon}, {returns, omega/(1 - alpha - beta)})));
d2log_likelihood_dalpha_dalpha = diff(log_likelihood, alpha, 2);
d2log_likelihood_dbeta_dbeta = diff(log_likelihood, beta, 2);
d2log_likelihood_dalpha_dbeta = diff(diff(log_likelihood, alpha),beta);
d2log_likelihood_dalpha_domega = diff(diff(log_likelihood, alpha),omega);
d2log_likelihood_dbeta_domega = diff(diff(log_likelihood, beta),omega);
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
estimated_parameters = [estimatedAlpha, estimatedBeta];
mle_FI=subs(Hessian_Matrix, [alpha, beta,omega], [estimatedAlpha, estimatedBeta,estimatedOmega]);
%cov0 = double(-inv(mle_FI));
cov0 = EstParamCov(2:3,2:3);

%% 1.Estimate the GARCH(1,1) model with a fixed constant (omega)
setOmega =0.07;
Mdl = garch('Constant',setOmega,'GARCHLags',1,'ARCHLags',1,'Offset',NaN);
[EstMdl,EstParamCov] = estimate(Mdl, returns);
disp('Estimated GARCH(1,1) Parameters:');
disp(EstMdl);
% Extract parameter estimates for alpha and beta
estimatedAlpha = EstMdl.ARCH{1};
estimatedBeta = EstMdl.GARCH{1};
estimated_parameters = [estimatedAlpha, estimatedBeta];
mle_FI=subs(Hessian_Matrix, [alpha, beta,omega], [estimatedAlpha, estimatedBeta,setOmega]);
%cov = double(-inv(mle_FI));
cov = EstParamCov(2:3,2:3);
%% 2.Estimate the GARCH(1,1) model with uncertained constant (omega)
Mdl2 = garch('GARCHLags',1,'ARCHLags',1,'Offset',NaN);
[EstMdl2,EstParamCov2] = estimate(Mdl2, returns);
estimatedAlpha2 = EstMdl2.ARCH{1};
estimatedBeta2 = EstMdl2.GARCH{1};
estimated_parameters2 = [estimatedAlpha2, estimatedBeta2];
estimatedconstant2 = EstMdl2.Constant;
mle_FI2=subs(Hessian_Matrix, [alpha, beta,omega], [estimatedAlpha2, estimatedBeta2,estimatedconstant2]);
%cov2 = double(-inv(mle_FI2));
cov2 = EstParamCov2(2:3,2:3);

%% adjusted cov for 1
%estimatedOmega=sigma_estimated*(1-estimatedAlpha-estimatedBeta);
c1=double(subs([d2log_likelihood_dalpha_domega;d2log_likelihood_dbeta_domega],[alpha, beta,omega], [estimatedAlpha, estimatedBeta, estimatedOmega]));
c=double(subs([d2log_likelihood_domega;d2log_likelihood_domega],[alpha, beta,omega], [estimatedAlpha, estimatedBeta, estimatedOmega]));
D_sym=cov0*[d2log_likelihood_dalpha_domega;d2log_likelihood_dbeta_domega];
D=double(subs(D_sym,[alpha, beta,omega], [estimatedAlpha, estimatedBeta,estimatedOmega]));
%v=variance_e_y_squared*(1-estimatedAlpha-estimatedBeta)^2;
%v=variance_e_y_squared*(1-estimated_nui_Alpha-estimated_nui_Beta)^2;
adjust=D*vv*transpose(D);
new_cov=cov0+D*vv*transpose(D);

%% adjusted cov for 2
%estimatedOmega2=sigma_estimated*(1-estimatedAlpha2-estimatedBeta2);
c2=double(subs([d2log_likelihood_dalpha_domega;d2log_likelihood_dbeta_domega],[alpha, beta,omega], [estimatedAlpha2, estimatedBeta2, estimatedOmega]));
D_sym2=cov0*[d2log_likelihood_dalpha_domega;d2log_likelihood_dbeta_domega];
%D_sym2=cov0*[d2log_likelihood_domega;d2log_likelihood_domega];
D2=double(subs(D_sym2,[alpha, beta,omega], [estimatedAlpha2, estimatedBeta2, estimatedOmega]));
%v2=variance_e_y_squared*(1-estimatedAlpha2-estimatedBeta2)^2;
adjust2=D2*vv*transpose(D2);
new_cov2=cov0+D2*vv*transpose(D2);

%% Plot the covariance ellipse for 1
critical_value = chi2inv(0.90, 2);
figure;
hold on
plot_ellipse(critical_value, cov, estimated_parameters,'r');
plot_ellipse(critical_value, cov0, estimated_parameters,'g');
plot_ellipse(critical_value, new_cov, estimated_parameters,'b');
scatter(estimated_parameters(1), estimated_parameters(2), 'k', 'filled');
scatter(ar, g, 'green', 'filled' )
axis equal;
xlabel('Alpha');
ylabel('Beta');
legend('unadjusted','true var', 'adjusted','estimated','true');
title({['90% Confidence Sigma Ellipse for Alpha and Beta:'];['omega=',num2str(setOmega),', true omega=',num2str(estimatedOmega)]});
grid on;
hold off

%% Plot the covariance ellipse for 2
figure;
hold on
plot_ellipse(critical_value, cov2, estimated_parameters2,'r');
plot_ellipse(critical_value, cov0, estimated_parameters2,'g');
plot_ellipse(critical_value, new_cov2, estimated_parameters2,'b')
scatter (estimated_parameters2(1),  estimated_parameters2(2), 'k', 'filled');
scatter(ar, g, 'green', 'filled' )
legend('unadjusted','true var','adjusted','estimated','true');
axis equal;
xlabel('Alpha');
ylabel('Beta');
title({['90% Confidence Sigma Ellipse for Alpha and Beta:'];['omega also estimated']});
grid on;
hold off
% Plot the sigma ellipse
