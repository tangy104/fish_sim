% Robotic Fish Simulation based on Dynamic Mathematical Model


% Main function to simulate over a 10-second period without input arguments
function main_simulation()
    % Define physical constants and parameters
    rho = 1000; % Density of water (kg/m^3)
    Cap = 1.479; % Added mass coefficient for the peduncle
    Cdp = 2.551; % Damping coefficient for the peduncle
    Cdf = 2.031; % Damping effect
    Caf = 0.134; % Added mass effect
    Lf = 0.05;  % Fin length
    Lp = 0.035; % Length of the peduncle (m)
    mf = 0.496 * 10^(-3);  % Total mass of the fin in kg
    alpha = 0.1; % Angle parameter (radians)
    mp = 0.5; % Mass of the robot body (kg)
    a = 0.033; % Fish robot dimensions (m)
    LCGp = 0.01; % Centre of gravity position (m)
    
    % Define constant theta values
    theta1 = 0.3516;
    theta2 = 0.0619;
    theta3 = 0.2502 * 10^-4;
    theta4 = 7.9023;
    theta5 = 8.7189;
    theta6 = -0.1914;
    theta7 = -0.0074;
    Ip = 0.01; % Moment of inertia (kg*m^2)

    % Define constants for the time-varying phi
    B = deg2rad(12);  % Amplitude of oscillation in radians
    f = 2.5;          % Frequency in Hz
    delta = 0;        % Phase offset
    phibar = 0;       % Baseline angle

    % Initial conditions for x (state vector)
    x0 = zeros(6,1);  % Initialize velocities and other states as zero

    % Simulation time parameters
    tspan = [0, 10];  % Simulate for 10 seconds

    % Solve the system of differential equations
    [t, x] = ode45(@(t, x) fish_dynamics(t, x, B, f, delta, phibar, rho, Cap, Lp, Lf, mp, mf, a, LCGp, theta1, theta2, theta3, theta4, theta5, theta6, theta7, Ip, alpha), tspan, x0);

    % Plotting results
    figure;
    subplot(2,1,1);
    plot(t, x(:, 1), 'r', t, x(:, 2), 'g', t, x(:, 3), 'b');
    title('Velocities over Time');
    xlabel('Time (s)');
    ylabel('Velocities (m/s and rad/s)');
    legend('u', 'v', 'r');

    subplot(2,1,2);
    plot(x(:, 1), x(:, 2));
    title('2D Trajectory');
    xlabel('u (m/s)');
    ylabel('v (m/s)');
end

% Call the main simulation
main_simulation();


% Function to calculate the time derivative of the state vector
function xdot = fish_dynamics(t, x, B, f, delta, phibar, rho, Cap, Cdf, Lp, Lf, mp, mf, a, LCGp, theta1, theta2, theta3, theta4, theta5, theta6, theta7, Ip, alpha)
    %syms xi 
    phi = B * sin(2 * pi * f * t + delta) + phibar;

    % Define `phione` as fundamental mode shape using Galerkin FEA
    phione = @(xi) sin(pi * xi / Lf)/sin(pi);

    % Define fin and related parameters
    df = @(xi) (2/5) * xi + 0.024;  % Thickness distribution function
    KM = @(xi) (1/12) * 0.215 * 10^(-1) * df(xi);
    Af = calculate_Af(df, Lf);  % Calculate fin area
    p = calculate_p(mf, df, Af); % Define p as a function of xi

    K = calculate_K(KM, phione, p, Lf);
    R = calculate_R(rho, rfdot, n, df, Cdf, phione, p, Lp);
    P = calculate_P(phione, p, df, rho, Caf, Lf);
    M = calculate_M(phione, p, df, rho, Caf, Lf);
    G = calculate_G(Lp, phione, p, df, rho, Caf, Lf);
    Fptilda = calculate_Fptilda(rho, rpdot, n, dp, Cdp, Lp);
    Mptilda = calculate_Mptilda(rho, rpdot, n, dp, Cdp, Lp);
    Z0 = calculate_Z0(rho, Cap, Lp);
    Z1 = calculate_Z1(rho, Cap, Lp);
    Z2 = calculate_Z2(rho, Cap, Lp);
    D1 = calculate_D1(phi, Z0, Z1, mp, a, LCGp);
    D2 = calculate_D2(phi, Z0);
    D3 = calculate_D3(phi, Z0);
    D4 = calculate_D4(Z1, mp, LCGp);
    D5 = theta4;
    D6 = theta5 * alpha;
    D7 = theta6 * alpha;  % Adjust if 'r / V' is calculated

    A = calculate_A_matrix(P, M, G, D1, D2, D3, D4, phi, phibar, theta1, theta2, theta3, Ip, Z2, mp, LCGp, a);
    B_mat = calculate_B_matrix(P, D2, D3, D4, D5, D6, D7, phi, phibar, mp, LCGp, a);
    C = calculate_C_vector(G, D4, phi, phibar, Ip, Z2, mp, LCGp, a);
    F = calculate_F_vector(K, R, phi, phibar, Fptilda, Mptilda, Lp, LCGp, mp);

    xdot = A \ (B_mat * x + F + C);
end

% Function to calculate the total fin area Af
function Af = calculate_Af(df, Lf)
    Af = integral(df, 0, Lf); % df should be an anonymous function of xi
end

% Function to calculate mass per unit length p
function p = calculate_p(mf, df, Af)
    p = @(xi) (mf * df(xi)) / Af; % p is now a function of xi
end

% Function to calculate Fptilda
function Fptilda = calculate_Fptilda(rho, rpdot, n, dp, Cdp, Lp)
    integrand = @(sigma) -(1/2) * rho * (rpdot * n) * abs(rpdot * n) * dp(sigma) * Cdp;
    Fptilda = integral(integrand, 0, Lp);
end

% Function to calculate Mptilda
function Mptilda = calculate_Mptilda(rho, rpdot, n, dp, Cdp, Lp)
    integrand = @(sigma) -(1/2) * rho * (rpdot * n) * abs(rpdot * n) * dp(sigma) * Cdp * sigma;
    Mptilda = integral(integrand, 0, Lp);
end

% Function to calculate R
function R = calculate_R(rho, rfdot, n, df, Cdf, phione, p, Lp)
    integrand = @(sigma) -(1/2) * rho * (rfdot * n) * abs(rfdot * n) * df(sigma) * Cdf * (phione / p);
    R = integral(integrand, 0, Lp);
end

% Function to calculate K
function K = calculate_K(KM, phione, p, Lf)
    % Define a symbolic variable for differentiation
    syms xi_sym

    % Define the integrand as a symbolic expression
    integrand = diff(KM(xi_sym) * diff(phione(xi_sym), xi_sym, 2), xi_sym, 2) * (phione(xi_sym) / p(xi_sym));
    
    % Convert integrand to a function handle
    integrand_func = matlabFunction(integrand, 'Vars', xi_sym);

    % Calculate K using numerical integration
    K = integral(@(xi) integrand_func(xi), 0, Lf, 'ArrayValued', true);
end

% Function to calculate P
function P = calculate_P(phione, p, df, rho, Caf, Lf)
    integrand = @(xi) phione * (1 + (pi * rho * Caf * df^2) / (4 * p));
    P = integral(integrand, 0, Lf);
end

% Function to calculate M
function M = calculate_M(phione, p, df, rho, Caf, Lf)
    integrand = @(xi) phione^2 * (1 + (pi * rho * Caf * df^2) / (4 * p));
    M = integral(integrand, 0, Lf);
end

% Function to calculate G
function G = calculate_G(Lp, phione, p, df, rho, Caf, Lf)
    integrand = @(xi) (Lp + xi) * phione * (1 + (pi * rho * Caf * df^2) / (4 * p));
    G = integral(integrand, 0, Lf);
end


% Functions to calculate Z0, Z1, Z2
function Z0 = calculate_Z0(rho, Cap, Lp)
    dp = @(sigma) 0.03 - (10/35) * sigma;
    Z0 = (1/4) * pi * rho * Cap * integral(@(sigma) dp(sigma).^2, 0, Lp);
end

function Z1 = calculate_Z1(rho, Cap, Lp)
    dp = @(sigma) 0.03 - (10/35) * sigma;
    Z1 = (1/4) * pi * rho * Cap * integral(@(sigma) dp(sigma).^2 .* sigma, 0, Lp);
end

function Z2 = calculate_Z2(rho, Cap, Lp)
    dp = @(sigma) 0.03 - (10/35) * sigma;
    Z2 = (1/4) * pi * rho * Cap * integral(@(sigma) dp(sigma).^2 .* sigma.^2, 0, Lp);
end

% Functions to calculate D1-D7
function D1 = calculate_D1(phi, Z0, Z1, mp, a, LCGp)
    D1 = Z0 * a * cos(phi)^2 + mp * a + mp * LCGp * cos(phi) + Z1 * cos(phi);
end

function D2 = calculate_D2(phi, Z0)
    D2 = Z0 * cos(phi)^2;
end

function D3 = calculate_D3(phi, Z0)
    D3 = Z0 * sin(phi) * cos(phi);
end

function D4 = calculate_D4(Z1, mp, LCGp)
    D4 = Z1 + mp * LCGp;
end

% Function to calculate the A matrix
function A = calculate_A_matrix(P, M, G, D1, D2, D3, D4, phi, phibar, theta1, theta2, theta3, Ip, Z2, mp, LCGp, a)
    A = [theta1 + D3 * sin(phibar), -(mp + D2) * sin(phibar), D1 * sin(phibar), 0, 0, 0;
        -D3 * cos(phibar), theta2 + (mp + D2) * cos(phibar), -D1 * cos(phibar), 0, 0, 0;
        a * D3 * cos(phibar) + D4 * sin(phi), -D4 * cos(phi) - a * (mp + D2) * cos(phibar), ...
        theta3 + Ip + Z2 + mp * LCGp * (LCGp + a * cos(phi)) + a * Z1 * cos(phi) + a * D1 * cos(phibar), 0, 0, 0;
        0, 0, 0, 1, 0, 0;
        0, 0, 0, 0, 1, 0;
        P * sin(phi), -P * cos(phi), P * a * cos(phi) + G, 0, 0, M];
end

% Function to calculate the B matrix
function B = calculate_B_matrix(P, D2, D3, D4, D5, D6, D7, phi, phibar, mp, LCGp, a)
    alpha = 1; % Adjust based on input
    B = [D6 * sin(alpha) - D5 * cos(alpha), D6 * sin(alpha) - D5 * cos(alpha), ...
        (mp * LCGp * sin(phi) - a * D3) * sin(phibar), 0, (mp + D2) * sin(phibar), D3 * sin(phibar);
        -D6 * cos(alpha) - D5 * sin(alpha), -D6 * cos(alpha) - D5 * sin(alpha), ...
        -(mp * LCGp * sin(phi) - a * D3) * cos(phibar), 0, -theta1 - (mp + D2) * cos(phibar), -D3 * cos(phibar);
        D7, D7, a * (mp * LCGp * sin(phi) - a * D3) * cos(phibar) - a * D4 * sin(phi), theta1 - theta2, ...
        D4 * cos(phi) + a * (mp + D2) * cos(phibar), D4 * sin(phi) + a * D3 * cos(phibar);
        0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0;
        0, 0, -a * P * sin(phi), 0, P * cos(phi), P * sin(phi)];
end

% Function to calculate the C vector
function C = calculate_C_vector(G, D4, phi, phibar, Ip, Z2, mp, LCGp, a)
    C = [-D4 * cos(phi) * sin(phibar); D4 * cos(phi) * cos(phibar); ...
        -(Z2 + Ip + mp * LCGp^2 + a * D4 * cos(phi) * cos(phibar)); 0; 0; -G];
end

% Function to calculate the F vector
function F = calculate_F_vector(K, R, phi, phibar, Fptilda, Mptilda, Lp, LCGp, mp)
    global tau Vf r phidot Mf niudot;
    F = [tau * cos(phibar) + (Fptilda * cos(phi) - Vf * cos(phi) + mp * LCGp * sin(phi) * (2 * r * phidot + phidot^2)) * sin(phibar);
        tau * sin(phibar) - (Fptilda * cos(phi) - Vf * cos(phi) + mp * LCGp * sin(phi) * (2 * r * phidot + phidot^2)) * cos(phibar);
        -a * tau * sin(phibar) + Mf - Vf * Lp + Mptilda + a * (Fptilda * cos(phi) - Vf * cos(phi) + mp * LCGp * sin(phi) * (2 * r * phidot + phidot^2)) * cos(phibar);
        r;
        qdot;
        (M * (niudot^2)-K) * q + R];
end

