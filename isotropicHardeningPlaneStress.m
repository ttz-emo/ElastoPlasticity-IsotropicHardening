%{
Copyright (c) 2026 Daniel Keller
Author: Daniel Keller
Version: 1.0
%}

function [sigma_n1, eps_p_n1, alpha_n1, eps_zz_n1, C_ep] = isotropicHardeningPlaneStress(eps_n, E_in, ET_in, nu_in, sigma_y_in, eps_p_n, alpha_n)   
    %% material parameter
    % compute the number of integer digits of the input value
    if E_in == 0
        n = 1;
    else
        n = floor(log10(abs(E_in))) + 1;
    end
    
    Gamma = 10^(n-1);       % stress scaling factor leads to higher numerical stability
    E = (1/Gamma)*E_in;     % young's modulus
    ET = (1/Gamma)*ET_in;   % elastoplastic tangent modulus
    nu = nu_in;             % poisson ratio
    sigma_y = (1/Gamma)*sigma_y_in; % yield stress
    K_bar = (E*ET)/(E-ET);  % isotropic hardening modulus
    mu = E/(2*(1 + nu));    % shear modulus (Lamé constants)
    
    % elasticity tensor
    C = (E/(1-nu^2)).*[1, nu, 0
                       nu, 1, 0
                       0 , 0, (1-nu)/2];
    
    % elastic compliance tensor
    C_inv = (1/E).*[1, -nu, 0
                    -nu, 1, 0
                    0  , 0, 2*(1+nu)];
    
    % projection tensor
    P = (1/3).*[2, -1, 0
                -1, 2, 0
                0,  0, 6];
    
    
    %% calc trial stress
    % trial stress
    sigma_trial = C * (eps_n - eps_p_n);
    
    % calculate relative stress (becomes only significant during kinematic hardening)
    xi_trial = sigma_trial;
    xi_trial_xx = xi_trial(1);
    xi_trial_yy = xi_trial(2);
    xi_trial_xy = xi_trial(3);
    
    
    %% yield condition
    f_trial = sqrt(xi_trial'*P*xi_trial) - sqrt(2/3)*(sigma_y + K_bar*alpha_n);
    
    
    %% elastic step
    if f_trial <= 1e-7
        sigma_n1 = sigma_trial.*Gamma;  % stress
        eps_p_n1 = eps_p_n;             % 2D plane stress plastic strain
        alpha_n1 = alpha_n;             % equivalent plastic strain
        eps_zz_n1 = (-nu/(E*Gamma))*(sigma_n1(1) + sigma_n1(2)) - (eps_p_n1(1) + eps_p_n1(2)); % normal strain in z-direction
        C_ep = C.*Gamma;                % elastoplastic tangent moduli
        return;
    end
    
    
    %% plastic step
    
    % eta represents the relative stress in the eigenraum (eta-space)
    eta_trial_11 = (xi_trial_xx + xi_trial_yy)/(sqrt(2));
    eta_trial_22 = (-xi_trial_xx + xi_trial_yy)/(sqrt(2));
    eta_trial_12 = xi_trial_xy;
    
    % inital parameter for newton iteration
    deltaGamma_k = 0;   % incremental plastic multiplier
    tolerance = 1e-10;  % termination criterion
    maxNbrIter = 50;    % maximum number of iteration
    
    % inital square of yield function
    f_square_k  = calc_f_square(deltaGamma_k, alpha_n, eta_trial_11, eta_trial_22, eta_trial_12, E, nu, mu, K_bar, sigma_y);
    
    % newton iteration
    for k = 1:maxNbrIter
        % derivation of square yield function
        df_square_k = calc_df_square(deltaGamma_k, alpha_n, eta_trial_11, eta_trial_22, eta_trial_12, E, nu, mu, K_bar, sigma_y);
        
        % newton step
        deltaGamma_k1 = deltaGamma_k - f_square_k / df_square_k;
    
        % square of yield function
        f_square_k1  = calc_f_square(deltaGamma_k1, alpha_n, eta_trial_11, eta_trial_22, eta_trial_12, E, nu, mu, K_bar, sigma_y);
        
        % termination criterion
        if abs(f_square_k1) < tolerance
            deltaGamma_k = deltaGamma_k1;
            break;
        end
        
        % update paramter
        deltaGamma_k = deltaGamma_k1;
        f_square_k = f_square_k1;
    end
    
    % verify that newton iteration has converged
    if(k >= maxNbrIter)
        error('Newton iteration isotropic hardening did not converge within the specified number of iterations.');
    end
    
    % verify deltaGamma >= 0
    if(deltaGamma_k < 0)
        error('The plastic multiplier is negative and therefore outside the correct definition.');
    end
    
    % set incremental plastic multiplier after newton iteration
    deltaGamma = deltaGamma_k;
    
    % compute modified elastic tangent moduli
    A1 = 2*E*deltaGamma + 3;
    B1 = 3*nu + E*deltaGamma;
    C1 = (1/2)*(E*deltaGamma - 3*nu + 3);
    D1 = E/(E^2*(deltaGamma + (2-nu)/(E))^2 - (2*nu-1)^2);
    Xi = D1 * [A1, B1, 0
               B1, A1, 0
               0 , 0 , C1];
    
    % calculate square of the norm of the relative stress in the eigenraum
    f_bar_sq = calc_f_bar_square(deltaGamma, eta_trial_11, eta_trial_22, eta_trial_12, E, nu, mu);
    
    % effective scalar stress (norm of the relative stress)
    f_bar = sqrt(f_bar_sq);
    
    % calculate relative stress (becomes only significant during kinematic hardening)
    xi_n1 = Xi * C_inv * xi_trial;
    
    % compute elastoplastic tangent moduli
    A2 = Xi*P*xi_n1;
    B2 = xi_n1'*P*Xi*P*xi_n1;
    Theta1 = 1;
    dK_dalpha = K_bar;
    Theta2 = 1 - (2/3)*dK_dalpha*deltaGamma;
    beta_bar_n1 = (2/3) * (Theta1/Theta2) * (dK_dalpha*Theta1) * xi_n1'*P*xi_n1;
    dsigma_deps = Xi - (A2*A2')/(B2 + beta_bar_n1);
    
    % update specific parameter
    sigma_n1 = xi_n1.*Gamma;                         % stress
    eps_p_n1 = eps_p_n + deltaGamma*P*xi_n1;         % 2D plane stress plastic strain
    alpha_n1 = alpha_n + sqrt(2/3)*deltaGamma*f_bar; % equivalent plastic strain
    eps_zz_n1 = (-nu/(E*Gamma))*(sigma_n1(1) + sigma_n1(2)) - (eps_p_n1(1) + eps_p_n1(2)); % normal strain in z-direction
    C_ep = dsigma_deps.*Gamma;                       % elastoplastic tangent moduli
end



%% calculate auxiliary variables in the eigenraum (eta-space)
function f_square = calc_f_square(deltaGamma, alpha_n, eta_trial_11, eta_trial_22, eta_trial_12, E, nu, mu, K_bar, sigma_y)
    % square of yield function
    f_square = (1/2) * calc_f_bar_square(deltaGamma, eta_trial_11, eta_trial_22, eta_trial_12, E, nu, mu) ...
                - calc_R_square(deltaGamma, alpha_n, eta_trial_11, eta_trial_22, eta_trial_12, E, nu, mu, K_bar, sigma_y);
end


function f_bar_square = calc_f_bar_square(deltaGamma, eta_trial_11, eta_trial_22, eta_trial_12, E, nu, mu)
    % define necessary constants    
    phi = E/(3*(1-nu));
    theta = 2*mu;

    % square of the norm of the relative stress in the eigenraum
    f_bar_square = (1/3) * (eta_trial_11^2)/((1+phi*deltaGamma)^2) + (eta_trial_22^2 + 2*eta_trial_12^2)/((1+theta*deltaGamma)^2);
end


function R_square = calc_R_square(deltaGamma, alpha_n, eta_trial_11, eta_trial_22, eta_trial_12, E, nu, mu, K_bar, sigma_y)
    % double squared radius of the yield surface (radius of yield surface is R = sqrt(2)*R(deltaGamma))
    R_square = (1/3) * (sigma_y + K_bar*(alpha_n + sqrt(2/3)*deltaGamma*abs(sqrt(calc_f_bar_square(deltaGamma, eta_trial_11, eta_trial_22, eta_trial_12, E, nu, mu)))))^2;
end


%% calculation of the derivative of auxiliary variables in the eigenraum (eta-space)
function df_square = calc_df_square(deltaGamma, alpha_n, eta_trial_11, eta_trial_22, eta_trial_12, E, nu, mu, K_bar, sigma_y)
    % derivation of square yield function
    df_square = (1/2) * calc_df_bar_square(deltaGamma, eta_trial_11, eta_trial_22, eta_trial_12, E, nu, mu) ...
                - calc_dR_square(deltaGamma, alpha_n, eta_trial_11, eta_trial_22, eta_trial_12, E, nu, mu, K_bar, sigma_y);
end


function df_bar_square = calc_df_bar_square(deltaGamma, eta_trial_11, eta_trial_22, eta_trial_12, E, nu, mu)
    % define necessary constants
    phi = E/(3*(1-nu));
    theta = 2*mu;

    % calculate derivation of f_bar_square
    df_bar_square = -((2*theta*(2*eta_trial_12^2 + eta_trial_22^2))/((deltaGamma*theta + 1)^3)) - ((2*eta_trial_11^2*phi)/(3*(deltaGamma*phi + 1)^3));    
end


function dR_square = calc_dR_square(deltaGamma, alpha_n, eta_trial_11, eta_trial_22, eta_trial_12, E, nu, mu, K_bar, sigma_y)
    % square of the norm of the relative stress in the eigenraum
    f_bar_sq = calc_f_bar_square(deltaGamma, eta_trial_11, eta_trial_22, eta_trial_12, E, nu, mu);

    % effective scalar stress (norm of the relative stress)
    f_bar = sqrt(f_bar_sq);

    % derivative of the squared norm of relative stress
    df_bar_sq = calc_df_bar_square(deltaGamma, eta_trial_11, eta_trial_22, eta_trial_12, E, nu, mu);
    
    % derivative of the effective scalar stress f_bar
    df_bar = df_bar_sq / (2 * f_bar);
    
    % current flow stress (radius) S = (sigma_y + K_bar * alpha_n1)
    S = sigma_y + K_bar * (alpha_n + sqrt(2/3) * deltaGamma * f_bar);
    
    % derivative of the flow stress term S
    dS = K_bar * sqrt(2/3) * (f_bar + deltaGamma * df_bar);
    
    % derivative of double squared radius of the yield surface (radius of yield surface is R = sqrt(2)*R(deltaGamma))
    dR_square = (2/3) * S * dS;    
end

