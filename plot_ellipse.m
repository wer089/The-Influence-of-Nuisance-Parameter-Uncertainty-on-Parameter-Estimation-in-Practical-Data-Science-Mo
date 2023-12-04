function plot_ellipse(critical_value, Fisher_Information_Matrix, parameter_estimates,color)
    % Determine the size and orientation of the ellipse using the Fisher Information Matrix
    % parameter_estimates: A 2-element vector with parameter estimates for alpha and beta
    
    % Calculate the standard errors for alpha and beta based on the Fisher Information Matrix
    standard_errors = sqrt(diag(inv(Fisher_Information_Matrix)));
    
    % Determine the lengths of the semi-major and semi-minor axes for the ellipse
    semi_major = sqrt(critical_value * Fisher_Information_Matrix(1, 1));
    semi_minor = sqrt(critical_value * Fisher_Information_Matrix(2, 2));
    
    % Compute the angle of rotation for the ellipse
    angle = atan2(Fisher_Information_Matrix(2, 1), Fisher_Information_Matrix(1, 1));
    
    % Create a range of angles for the ellipse
    theta = linspace(0, 2 * pi, 100);
    
    % Parametric equations for the ellipse
    x = semi_major * cos(theta);
    y = semi_minor * sin(theta);
    
    % Rotate the ellipse by the calculated angle
    ellipse_x = x * cos(angle) - y * sin(angle);
    ellipse_y = x * sin(angle) + y * cos(angle);
    
    % Translate the ellipse to the parameter estimates
    ellipse_x = ellipse_x + parameter_estimates(1);
    ellipse_y = ellipse_y + parameter_estimates(2);
    
    % Plot the sigma ellipse
    plot(ellipse_x, ellipse_y, color);
    %scatter(parameter_estimates(1), parameter_estimates(2), 'k', 'filled');  % Parameter estimates
    
end
