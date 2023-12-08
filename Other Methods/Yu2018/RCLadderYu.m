function [sys_est, RC_est] = RCLadderYu(R_scaled, C_scaled, T, is_prescale)

    %% Change the parameterization:
    R_scaler = 10^-floor(mean(log10(R_scaled)));
    C_scaler = 10^-floor(mean(log10(C_scaled)));
    R_scaled = R_scaled*R_scaler;
    C_scaled = C_scaled*C_scaler;
    RC2theta = @(R, C) 1./(repelem(R, 2).*[repelem(C, [1 2*ones(1, length(C)-1)]); 1]);
    
    fn_A = @(theta) -diag(theta(1:2:end)) - diag([0; theta(2:2:end)]) + diag(theta(1:2:end-2), 1) + diag(theta(2:2:end), -1);
    fn_B = @(theta) [zeros((length(theta)-1)/2, 1); theta(end)];
    fn_C = @(theta) [zeros(1, (length(theta)-1)/2) 1/R_scaled(end)];
    fn_D = @(theta) -1/R_scaled(end);

    theta = RC2theta(R_scaled, C_scaled);
    theta = theta(1:end-1); % Last element is known from D.
    sys = ss(fn_A(theta), fn_B(theta), fn_C(theta), fn_D(theta), 1);
    sys = ss2ss(sys, T);
    if is_prescale
        sys = prescale(sys);
    end
    
    %% Run the optimization
    [sys_est, theta_est] = Yu2018(sys, fn_A, fn_B, fn_C, length(theta), 'verbose', 0);
    sys_est = ss(sys_est.A, sys_est.B, sys_est.C, -1/R_scaled(end));
    [R_est_scaled, C_est_scaled] = theta2RC(theta_est, 1/R_scaled(end));
    RC_est = [R_est_scaled/R_scaler; C_est_scaled/C_scaler];
    
    function [R, C] = theta2RC(theta, invRend)
        theta = [theta; invRend];
        theta_inv = 1./theta;
        
        R = zeros(length(theta), 1);
        R(1:2:end) = theta_inv(end:-2:1);
        R(2:2:end) = theta(end-1:-2:1);
        R = cumprod(R);
        R = R(end-1:-2:1);
        
        C = zeros(length(theta), 1);
        C(1:2:end) = theta(end:-2:1);
        C(2:2:end) = theta_inv(end-1:-2:1);
        C = cumprod(C);
        C = C(end:-2:1);
    end
    
end