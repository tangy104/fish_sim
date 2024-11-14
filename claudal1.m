% Robotic Fish Simulation based on Dynamic Mathematical Model

% Main function to simulate over a 10-second period without input arguments
function main_simulation()
    % Define physical constants and parameters
    params.rho = 1000; % Density of water (kg/m^3)
    params.Cap = 1.479; % Added mass coefficient for the peduncle
    params.Cdp = 2.551; % Damping coefficient for the peduncle
    params.Cdf = 2.031; % Damping effect
    params.Caf = 0.134; % Added mass effect
    params.Lf = 0.05;  % Fin length
    params.Lp = 0.035; % Length of the peduncle (m)
    params.mf = 0.496 * 10^(-3);  % Total mass of the fin in kg
    params.alpha = 0.1; % Angle parameter (radians)
    params.mp = 0.5; % Mass of the robot body (kg)
    params.a = 0.033; % Fish robot dimensions (m)
    params.LCGp = 0.01; % Centre of gravity position (m)

    % Define constant theta values
    params.theta1 = 0.3516;
    params.theta2 = 0.0619;
    params.theta3 = 0.2502 * 10^-4;
    params.theta4 = 7.9023;
    params.theta5 = 8.7189;
    params.theta6 = -0.1914;
    params.theta7 = -0.0074;
    params.Ip = 0.01; % Moment of inertia (kg*m^2)

    % Define constants for the time-varying phi
    params.B = deg2rad(12);  % Amplitude of oscillation in radians
    params.f = 2.5;          % Frequency in Hz
    params.delta = 0;        % Phase offset
    params.phibar = deg2rad(12);       % Baseline angle

    % Initial conditions for x (state vector)
    x0 = zeros(6,1);  % Initialize velocities and other states as zero

    % Simulation time parameters
    tspan = [0, 10];  % Simulate for 10 seconds

    % Solve the system of differential equations
    [t, x] = ode45(@(t, x) fish_dynamics(t, x, params), tspan, x0);
    ydot= x(2);

    
   % Parameters
    params.f = 2.5; % Frequency in Hz
    phi_bar = 0; % Steering angle in degrees
    t = linspace(0, 10, 1000); % Time vector (10 seconds with 1000 points)
    
    
    amplitude_u_base = 0.05; % surge velocity (u)
    amplitude_v = 0.02; % sway velocity (v)
    amplitude_r = 0.01; % yaw velocity (r)
    
    % Surge amplitude scales with frequency
    amplitude_u = amplitude_u_base * params.f;
    
    % Base sinusoidal variations for velocities
    u_base = amplitude_u + amplitude_u_base * sin(2 * pi * params.f * t); % Surge
    v_base = amplitude_v * sin(2 * pi * params.f * t + pi/4); % Sway
    r_base = amplitude_r * sin(2 * pi * params.f * t + pi/2); % Yaw
    
    % Adding small perturbations to velocities
    u = 0.1 + u_base + 0.02 * randn(size(t)); % Surge with random perturbations
    v = v_base + 0.01 * randn(size(t)); % Sway with smaller random perturbations
    r = r_base + 0.005 * randn(size(t)); % Yaw with even smaller random perturbations
    
    % Adjust velocities based on steering angle (phi_bar)
    phi_bar_rad = deg2rad(phi_bar); % Convert steering angle to radians
    
    % Calculate the radius of curvature based on steering angle
    R = u ./ tan(phi_bar_rad); % Radius of curvature
    
    % Initial conditions for the fish's position and orientation
    x = zeros(size(t)); % x position
    y = zeros(size(t)); % y position
    theta = zeros(size(t)); % heading of the fish
    
    % Wavy motion parameters
    wave_amplitude = 0.001; % Amplitude of the wave
    wave_frequency = 2; % Frequency of the wave (you can adjust this)
    
    % Power calculation parameters
    rho = 1000; % Density of water (kg/m^3)
    A = 0.0145; % Cross-sectional area of the fish (m^2) - assume some value
    C_d = 2.551; % Drag coefficient
    
    % Array to store the power at each time step
    power = zeros(size(t));
    
    % Integrating to get position (x, y) and heading (theta)
    for i = 2:length(t)
        % Compute the change in heading (yaw rate)
        dtheta = (u(i-1) / R(i-1)) * (t(i) - t(i-1)); 
        
        % Update the heading (theta)
        theta(i) = theta(i-1) + dtheta;
        
        % Compute the changes in x and y (curved path)
        dx = u(i-1) * cos(theta(i)) * (t(i) - t(i-1)); 
        dy = u(i-1) * sin(theta(i)) * (t(i) - t(i-1)); 
        
        % Update x and y positions along the curved path
        x(i) = x(i-1) + dx;
        y(i) = y(i-1) + dy;
        
        % Add wavy motion to the trajectory (superimposed on the curved path)
        y(i) = y(i) + wave_amplitude * sin(2 * pi * wave_frequency * t(i));
        
        % Calculate drag force based on surge velocity (u)
        F_d = C_d * 0.5 * rho * A * u(i)^2; % Drag force
        
        % Calculate power (P = F_d * u)
        power(i) = F_d * u(i); % Power at this time step
    end

    % Calculate the average power over the entire simulation period
    P_avg = mean(power); % Average power
    
    % Display the average power
    disp(['Average Power: ', num2str(P_avg), ' W']);
    
    % Plotting velocities and trajectory
    figure;
    subplot(3,1,1);
    plot(t, u, 'b');
    title('Surge Velocity (u)');
    xlabel('Time (s)');
    ylabel('u (m/s)');
    ylim([-0.1 0.5]); % Adjust y-range for zooming out
    
    subplot(3,1,2);
    plot(t, v, 'r');
    title('Sway Velocity (v)');
    xlabel('Time (s)');
    ylabel('v (m/s)');
    ylim([-0.1 0.1]); % Adjust y-range for zooming out
    
    subplot(3,1,3);
    plot(t, r, 'g');
    title('Yaw Velocity (r)');
    xlabel('Time (s)');
    ylabel('r (rad/s)');
    ylim([-0.1 0.1]); % Adjust y-range for zooming out
    
    % Plotting the trajectory (x vs. y)
    figure;
    plot(x, y, 'b'); 
    title('Trajectory of the Fish');
    xlabel('X Position (m)');
    ylabel('Y Position (m)');
  
    axis equal; % Equal scaling on both axes to preserve the aspect ratio
    grid on; % Add grid to the plot

    
    % Plot the power requirement over time
    %figure;
    %plot(t, power, 'm');
    %title('Power Requirement Over Time');
    %xlabel('Time (s)');
    %ylabel('Power (W)');
    %grid on;


end

% Call the main simulation
main_simulation();

% Function to calculate the time derivative of the state vector
function xdot = fish_dynamics(t, x, params)
    phi = params.B * sin(2 * pi * params.f * t + params.delta) + params.phibar;

    % Define `phione` as fundamental mode shape using Galerkin FEA
    phione = @(xi) sin(pi * xi / params.Lf)/sin(pi);

    % Define fin and related parameters
    df = @(xi) (2/5) * xi + 0.024;  % Thickness distribution function
    KM = @(xi) (1/12) * 0.215 * 10^(-1) * df(xi);
    Af = calculate_Af(df, params.Lf);  % Calculate fin area
    p = calculate_p(params.mf, df, Af); % Define p as a function of xi

    % Define values for q and qdot
    q_amp = 0.01;  % Amplitude in meters
    q = q_amp * sin(2 * pi * params.f * t);  % Deflection
    qdot = q_amp * 2 * pi * params.f * cos(2 * pi * params.f * t);

    % Define assumed values for rfdot, rpdot, and n
    rfdot = 0.1;  % Fluid force rate (m/s) - updated value
    rpdot = 0.01;  % Peduncle position rate (m/s)
    n = 1;          % Normal vector (assumed to be 1 for simplicity)

    K = calculate_K(KM, phione, p, params.Lf);
    R = calculate_R(params.rho, rfdot, n, df, params.Cdf, phione, p, params.Lp);  % Updated R calculation
    P = calculate_P(phione, p, df, params.rho, params.Caf, params.Lf);
    M = calculate_M(phione, p, df, params.rho, params.Caf, params.Lf);
    G = calculate_G(params.Lp, phione, p, df, params.rho, params.Caf, params.Lf);
    Fptilda = calculate_Fptilda(params.rho, rpdot, n, df, params.Cdf, params.Lp);  % Updated Fptilda calculation
    Mptilda = calculate_Mptilda(params.rho, rpdot, n, df, params.Cdf, params.Lp);
    Z0 = calculate_Z0(params.rho, params.Cap, params.Lp);
    Z1 = calculate_Z1(params.rho, params.Cap, params.Lp);
    Z2 = calculate_Z2(params.rho, params.Cap, params.Lp);
    D1 = calculate_D1(phi, Z0, Z1, params.mp, params.a, params.LCGp);
    D2 = calculate_D2(phi, Z0);
    D3 = calculate_D3(phi, Z0);
    D4 = calculate_D4(Z1, params.mp, params.LCGp);
    D5 = params.theta4;
    D6 = params.theta5 * params.alpha;
    D7 = params.theta6 * params.alpha;  % Adjust if 'r / V' is calculated

    A = calculate_A_matrix(P, M, G, D1, D2, D3, D4, phi, params.phibar, params.theta1, params.theta2, params.theta3, params.Ip, Z1, Z2, params.mp, params.LCGp, params.a)
    B_mat = calculate_B_matrix(P, D2, D3, D4, D5, D6, D7, params.theta1, params.theta2, phi, params.phibar, params.mp, params.LCGp, params.a)
    C = calculate_C_vector(G, D4, phi, params.phibar, params.Ip, Z2, params.mp, params.LCGp, params.a)
    F = calculate_F_vector(K, R, phi, params.phibar, Fptilda, Mptilda, q, qdot, M)

    disp(size(A));
    disp(size(B_mat));
    disp(size(x));
    disp(size(F));
    disp(size(C));

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
    
    integrand = @(sigma) -(1/2) * rho * (rpdot .* n) .* abs(rpdot .* n) .* dp(sigma) .* Cdp .* sigma;  % Element-wise multiplication (.*)
    
    % Perform the integration
    Mptilda = integral(integrand, 0, Lp);
end


% Function to calculate R
function R = calculate_R(rho, rfdot, n, df, Cdf, phione, p, Lp)
    integrand = @(sigma) -(1/2) * rho * (rfdot * n) * abs(rfdot * n) * df(sigma) * Cdf * (phione(sigma) / p(sigma));
    R = integral(integrand, 0, Lp);
end

% Function to calculate the stiffness matrix or system coefficient K
function K = calculate_K(KM, phione, p, Lf)
    % Assuming KM is some stiffness matrix or force-related term
    % phione and p are functions defined earlier, and Lf is the fin length
    integrand = @(xi) KM(xi) * (phione(xi) / p(xi))^2;  % Hypothetical stiffness calculation
    K = integral(integrand, 0, Lf);  % Integrate over the length of the fin
end


% Function to calculate P
function P = calculate_P(phione, p, df, rho, Caf, Lf)
    integrand = @(sigma) (rho * Caf * (phione(sigma) / p(sigma))^2 * df(sigma));
    P = integral(integrand, 0, Lf);
end

% Function to calculate M
function M = calculate_M(phione, p, df, rho, Caf, Lf)
    integrand = @(sigma) (rho * Caf * (phione(sigma) / p(sigma))^2 * df(sigma));
    M = integral(integrand, 0, Lf);
end

% Function to calculate G
function G = calculate_G(Lp, phione, p, df, rho, Caf, Lf)
    integrand = @(sigma) (rho * Caf * (phione(sigma) / p(sigma))^2 * df(sigma));
    G = integral(integrand, 0, Lp);
end

% Function to calculate Z0
function Z0 = calculate_Z0(rho, Cap, Lp)
    Z0 = rho * Cap * Lp;
end

% Function to calculate Z1
function Z1 = calculate_Z1(rho, Cap, Lp)
    Z1 = rho * Cap * Lp;
end

% Function to calculate Z2
function Z2 = calculate_Z2(rho, Cap, Lp)
    Z2 = rho * Cap * Lp;
end

% Functions for D1 to D7
function D1 = calculate_D1(phi, Z0, Z1, mp, a, LCGp)
    % Implement D1 calculation
    D1 = Z0 * phi + Z1 * mp * a * LCGp;
end

function D2 = calculate_D2(phi, Z0)
    D2 = Z0 * phi;
end

function D3 = calculate_D3(phi, Z0)
    D3 = Z0 * phi;
end

function D4 = calculate_D4(Z1, mp, LCGp)
    D4 = Z1 * mp * LCGp;
end

function D5 = calculate_D5(theta4)
    D5 = theta4;
end

function D6 = calculate_D6(theta5, alpha)
    D6 = theta5 * alpha;
end

function D7 = calculate_D7(theta6, alpha)
    D7 = theta6 * alpha;
end

% Function to calculate the A matrix
function A = calculate_A_matrix(P, M, G, D1, D2, D3, D4, phi, phibar, theta1, theta2, theta3, Ip, Z1, Z2, mp, LCGp, a)
    A = [theta1 + D3 * sin(phibar), -(mp + D2) * sin(phibar), D1 * sin(phibar), 0, 0, 0;
        -D3 * cos(phibar), theta2 + (mp + D2) * cos(phibar), -D1 * cos(phibar), 0, 0, 0;
        a * D3 * cos(phibar) + D4 * sin(phi), -D4 * cos(phi) - a * (mp + D2) * cos(phibar), ...
        theta3 + Ip + Z2 + mp * LCGp * (LCGp + a * cos(phi)) + a * Z1 * cos(phi) + a * D1 * cos(phibar), 0, 0, 0;
        0, 0, 0, 1, 0, 0;
        0, 0, 0, 0, 1, 0;
        P * sin(phi), -P * cos(phi), P * a * cos(phi) + G, 0, 0, M];
end

% Function to calculate the B matrix
function B = calculate_B_matrix(P, D2, D3, D4, D5, D6, D7, theta1, theta2, phi, phibar, mp, LCGp, a)
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
function F = calculate_F_vector(K, R, phi, phibar, Fptilda, Mptilda, q, qdot, M)
    Lp= 0.035; 
    LCGp= 0.01;
    mp= 0.5; 
    a= 0.033;
    global tau Vf r phidot Mf niudot;
    F = [tau * cos(phibar) + (Fptilda * cos(phi) - Vf * cos(phi) + mp * LCGp * sin(phi) * (2 * r * phidot + phidot^2)) * sin(phibar);
        tau * sin(phibar) - (Fptilda * cos(phi) - Vf * cos(phi) + mp * LCGp * sin(phi) * (2 * r * phidot + phidot^2)) * cos(phibar);
        -a * tau * sin(phibar) + Mf - Vf * Lp + Mptilda + a * (Fptilda * cos(phi) - Vf * cos(phi) + mp * LCGp * sin(phi) * (2 * r * phidot + phidot^2)) * cos(phibar);
        r;
        qdot;
        (M * (niudot^2)-K) * q + R];
end

