%{
Copyright (c) 2026 Daniel Keller
Author: Daniel Keller
Version: 1.0

This code performs a linear isotropic plastic hardening calculation for the two-dimensional state plane stress.

main literature:
[1] J. C. Simo und T. J. R. Hughes, Hrsg., „Computational Inelasticity“, in Computational Inelasticity, New York, NY: Springer, 1998, S. 1–70. doi: 10.1007/0-387-22763-6_1.
[2] N.-H. Kim, Introduction to Nonlinear Finite Element Analysis. New York, NY: Springer US, 2015. doi: 10.1007/978-1-4419-1746-1.
[3] W. Rust, Nichtlineare Finite-Elemente-Berechnungen: Kontakt, Geometrie, Material. Wiesbaden: Vieweg+Teubner, 2011. doi: 10.1007/978-3-8348-8148-9.
[4] M. Wächter, C. Müller, und A. Esderts, Angewandter Festigkeitsnachweis nach FKM-Richtlinie: Kurz und bündig. Wiesbaden: Springer Fachmedien Wiesbaden, 2021. doi: 10.1007/978-3-658-32857-3.
%}
clear;
close all;

% add all folders to path
currentFolderPath = pwd();
addpath(genpath(currentFolderPath));

% save figures in current folder
saveFig = false;


%% simulation parameter
% simulated strain excitation
% normal_x: strain in x-direction
% normal_y: strain in y-direction
% shear: strain in xy-direction
% combined: strain in x, y and xy-direction
forceCase = "normal_x";   

epsMax = 1e-2;  % maximum compoment strain in Voigt notation for the selected case
epsStep = 1e-5; % step size for compoment strain in Voigt notation for the selected case


%% material parameter
E = 200e9;       % young's modulus
ET = 1.45e9;     % elastoplastic tangent modulus
nu = 0.3;        % poisson ratio
sigma_y = 250e6; % yield stress


%% simulated strain excitation
switch forceCase
    case "normal_x"
        plotText = 'strain excitation in x-direction';
        eps = [0:epsStep:epsMax; zeros(2, numel(0:epsStep:epsMax))];
    case "normal_y"
        plotText = 'strain excitation in y-direction';
        eps = [zeros(1, numel(0:epsStep:epsMax)); 0:epsStep:epsMax; zeros(1, numel(0:epsStep:epsMax))];
    case "shear"
        plotText = 'strain excitation in xy-direction';
        eps = [zeros(2, numel(0:epsStep:epsMax)); 0:epsStep:epsMax];
    case "combined"
        plotText = 'combined strain excitation in x, y and xy-direction';
        eps = [0:epsStep:epsMax; 0:epsStep:epsMax; 0:epsStep:epsMax];
end


%% FEM calculation (example for one gauss point)
% values for plastic calculation
eps_p_n = [0, 0, 0]'; % intial plane stress plastic strain
alpha_n = 0;          % inital equivalent plastic strain

% values to storage
eps_total_VM = zeros(1,length(eps(1,:)));
eps_p_VM = zeros(1,length(eps(1,:)));
eps_e_VM = zeros(1,length(eps(1,:)));
sigma_VM = zeros(1,length(eps(1,:)));

for i = 1:length(eps(1,:))
    % set current strain tensor
    eps_n = eps(:,i);

    % calculation performed for each Gauss point
    [sigma_n1, eps_p_n1, alpha_n1, eps_zz_n1, C_ep] = isotropicHardeningPlaneStress(eps_n, E, ET, nu, sigma_y, eps_p_n, alpha_n);
    
    sigma_n = sigma_n1;   % stress
    eps_p_n = eps_p_n1;   % plane stress plastic strain
    alpha_n = alpha_n1;   % equivalent plastic strain
    eps_zz_n = eps_zz_n1; % normal strain in z-direction
        
    % plastic equivalent von Mises strain (Attention: In J2 plasticity the effective plastic Poisson's ratio in nu_p = 0.5 because plastic flow is assumed to be isochoric, see [1, p. 95].)
    eps_p_VM_n = alpha_n;
    
    % elastic equivalent von Mises strain
    eps_e_n = eps_n - eps_p_n;
    eps_xx = eps_e_n(1); eps_yy = eps_e_n(2); eps_zz = -(nu)/(1-nu)*(eps_e_n(1) + eps_e_n(2)); gamma_xy = eps_e_n(3);
    eps_e_VM_n = vonMisesStrainPlaneStress(eps_xx, eps_yy, eps_zz, gamma_xy, nu);

    % total equivalent von Mises strain (This is only an approximation. From a mathematical perspective, the plastic and elastic von Mises strains should not be added together directly.)
    eps_total_VM_n = eps_e_VM_n + eps_p_VM_n;
    
    % equivalent von Mises stress
    sigma_xx = sigma_n(1); sigma_yy = sigma_n(2); sigma_xy = sigma_n(3);
    sigma_VM_n = vonMisesStressPlaneStress(sigma_xx, sigma_yy, sigma_xy);

    % storage values
    eps_total_VM(i) = eps_total_VM_n;
    eps_p_VM(i) = eps_p_VM_n;
    eps_e_VM(i) = eps_e_VM_n;
    sigma_VM(i) = sigma_VM_n;
end


%% plot
% get file path
filePath = fileparts(mfilename('fullpath'));

% expanded elastic domain
epsNewElasticityLine = [0, eps_e_VM_n] + (eps_total_VM_n - eps_e_VM_n);
sigmaNewElasticityLine = [0, eps_e_VM_n*E];

% plot stress-strain curve
namePlot = [plotText, ' stress strain curve'];
fig = figure('Name', namePlot, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
hold on
plot(eps_total_VM*100, sigma_VM*1e-6, '-', 'LineWidth', 4, 'DisplayName', 'stress-strain curve');
plot(epsNewElasticityLine*100, sigmaNewElasticityLine*1e-6, '-', 'LineWidth', 4, 'DisplayName', 'expanded elastic domain');
hold off
xlabel('total von Mises strain $\varepsilon_\mathrm{total,VM}$ in \%', 'Interpreter', 'latex');
ylabel('von Mises stress $\sigma_\mathrm{VM}$ in MPa', 'Interpreter', 'latex');
title(['linear isotropic hardening plane stress: ', plotText], 'Interpreter', 'latex');
grid on;
legend('Interpreter', 'latex', 'FontSize', 15, 'NumColumns', 1);
set(gca, 'FontSize', 20);
if(saveFig)
    exportgraphics(fig, [filePath, '\', namePlot, '.png'], 'Resolution', 300);
    savefig(fig, [filePath, '\', namePlot, '.fig'])
end

% plot
namePlot = [plotText, ' strain curves'];
fig = figure('Name', namePlot, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
hold on
plot(eps_total_VM, eps_total_VM, '-', 'LineWidth', 4, 'DisplayName', 'equivalent total strain');
plot(eps_total_VM, eps_e_VM, '-', 'LineWidth', 4, 'DisplayName', 'equivalent elastic strain');
plot(eps_total_VM, eps_p_VM, '-', 'LineWidth', 4, 'DisplayName', 'equivalent plastic strain');
hold off
xlabel('total von Mises strain $\varepsilon_\mathrm{total,VM}$ in \%', 'Interpreter', 'latex');
ylabel('von Mises strain $\varepsilon_\mathrm{VM}$ in \%', 'Interpreter', 'latex');
title(['linear isotropic hardening plane stress: ', plotText], 'Interpreter', 'latex');
grid on;
legend('Interpreter', 'latex', 'FontSize', 15, 'NumColumns', 1);
set(gca, 'FontSize', 20);
if(saveFig)
    exportgraphics(fig, [filePath, '\', namePlot, '.png'], 'Resolution', 300);
    savefig(fig, [filePath, '\', namePlot, '.fig'])
end


%% additional functions
% equivalent von Mises strain for plane stress
function eps_VM = vonMisesStrainPlaneStress(eps_xx, eps_yy, eps_zz, gamma_xy, nu)
    eps_VM = (1/(1+nu)) * sqrt((1/2) * ((eps_xx-eps_yy)^2 + (eps_yy-eps_zz)^2 + (eps_zz-eps_xx)^2 + (3/2)*gamma_xy^2));
end

% equivalent von Mises stress for plane stress
function sigma_VM = vonMisesStressPlaneStress(sigma_xx, sigma_yy, sigma_xy)
    sigma_VM = sqrt(sigma_xx^2 - sigma_xx*sigma_yy + sigma_yy^2 + 3*sigma_xy^2);
end

