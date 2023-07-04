%IDMONTECARLO 
%   Script to run the Monte Carlo simulations for Section XXX of [1]
%
%   See also: RCLADDERN, RCLADDERSTRUCTUREDID, RCLADDERHWANG,
%             RCLADDER2THETA
%       [1] Hannigan, B. C., Menon, C. (Draft). Fast, analytical method for 
%           structured identification of SISO RC-ladder-type systems.
%       [2] Hwang, C., Nakano, B., & Asahi, T. (1984). Transformation of 
%           state-space model to Cauer I and II CFE canonical forms. 
%           International Journal of Systems Science, 15(7), 797–804. 
%           https://doi.org/10.1080/00207728408926600
%
%   $Author: BH$    $Date: 2023-06-13$  $Revision: 1$
%
%   REVISION 1: 2023-06-28
%       Merged in code for simulation comparing the algorithm with IDGREY 
%       when identifying an RC-ladder system from frequency-domain data.
%
%   ©2023 ETH Zurich, Brett Hannigan; D-HEST; Biomedical and Mobile Health Technology (BMHT) Lab; Carlo Menon

fprintf('Script to reproduce Monte Carlo simulations from:\n')
fprintf('>\t Hannigan, B. C., Menon, C. (Draft). Fast, analytical method for structured identification of SISO RC-ladder-type systems.\n\n')

% True system parameters:
RANGE_R = [1 10];
RANGE_C = [1 10];
SCALING_R = 1e4; % x10 kOhm
SCALING_C = 1e-11; % pF

% Options for Section XXX
IS_PRESCALE = true; % Use MATLAB PRESCALE() function on randomly transformed systems.

% Options for Section YYY
freq_range = [1e3 1e7];
n_freq = 128; % Number of test points for frequency domain identification.
w_test = 2*pi*logspace(log10(freq_range(1)), log10(freq_range(2)), n_freq);

% Monte Carlo options:
order = 3:20;
n_trials = 100;
rng(1, 'twister') % Seed RNG for reproducibility.
s = struct('order', num2cell(repelem(order', n_trials))); % Results structure.

fprintf('Starting simulation of %i systems\n', length(order)*n_trials)
fprintf('Parameters:\n')
fprintf('\t R: U~(%.2e, %.2e) ohm', RANGE_R(1)*SCALING_R, RANGE_R(2)*SCALING_R)
fprintf('\t C: U~(%.2e, %.2e) F\n', RANGE_C(1)*SCALING_C, RANGE_C(2)*SCALING_C)
fprintf('\t Order: [%i, ..., %i] \t\t\t Trials: %i\n', order(1), order(end), n_trials);

for i_n=1:length(order)
    for i_trial=1:n_trials    
        %% Generate and Simulate Random True State-Space System.
        index = (i_n - 1)*n_trials + i_trial;
        fprintf('\n Running simulation %i/%i', index, length(order)*n_trials)
        n = order(i_n);
        R = (RANGE_R(1) + rand(n, 1)*diff(RANGE_R))*SCALING_R;
        C = (RANGE_C(1) + rand(n, 1)*diff(RANGE_C))*SCALING_C;
        s(index).sys_true = RCLadderN(R, C, 'ascending', false);
        s(index).R_true = R;
        s(index).C_true = C;
        
        %% Generate Random State Transformation and Transform System
        T = randn(n);
        s(index).T = T;
        if IS_PRESCALE
            s(index).sys_t = prescale(ss2ss(s(index).sys_true, T));
        else
            s(index).sys_t = ss2ss(s(index).sys_true, T);
        end
       
        %% Run Reconstruction Algorithm
        % The main algorithm from [1].
        try
            tic
            [sys_est, T_est] = RCLadderStructuredID(s(index).sys_t, 1);
            s(index).duration_est = toc;
            s(index).sys_est = sys_est;
            s(index).A_est_dist = MatrixDistance(s(index).sys_true.A, s(index).sys_est.A);
            s(index).B_est_dist = MatrixDistance(s(index).sys_true.B, s(index).sys_est.B);
            s(index).C_est_dist = MatrixDistance(s(index).sys_true.C, s(index).sys_est.C);
            s(index).D_est_dist = MatrixDistance(s(index).sys_true.D, s(index).sys_est.D);
        catch
            s(index).duration_est = nan;
            s(index).sys_est = [];
            s(index).A_est_dist = nan;
            s(index).B_est_dist = nan;
            s(index).C_est_dist = nan;
            s(index).D_est_dist = nan;
        end
        
        %% Run Routh Array Procedure
        % The procedure from [2] generates the structured state-space
        % matrices from the Routh array coefficients and Markov parameters
        % of a Cauer type I system.
        try
            tic
            [sys_hwang, ~] = RCLadderHwang(s(index).sys_t);
            s(index).duration_hwang = toc;
            s(index).sys_hwang = sys_hwang;
            s(index).A_hwang_dist = MatrixDistance(s(index).sys_true.A, s(index).sys_hwang.A);
            s(index).B_hwang_dist = MatrixDistance(s(index).sys_true.B, s(index).sys_hwang.B);
            s(index).C_hwang_dist = MatrixDistance(s(index).sys_true.C, s(index).sys_hwang.C);
            s(index).D_hwang_dist = MatrixDistance(s(index).sys_true.D, s(index).sys_hwang.D);
        catch
            s(index).duration_hwang = nan;
            s(index).sys_hwang = [];
            s(index).A_hwang_dist = nan;
            s(index).B_hwang_dist = nan;
            s(index).C_hwang_dist = nan;
            s(index).D_hwang_dist = nan;
        end

        %% Collect Frequency-Domain Data
        % Exact frequency response data (no added noise).
        [mag, phase, ~] = bode(s(index).sys_true, w_test);
        response = squeeze(mag).*exp(1i*deg2rad(squeeze(phase)));
        % Create IDFRD object out of frequency response data.
        freq_id = idfrd(response, w_test, 0);
        
        %% System Identification using IDGREY
        % Generate random initial guess within the range of the parameters.
        initial_theta = [RANGE_R(1) + rand(n, 1)*diff(RANGE_R); RANGE_C(1) + rand(n, 1)*diff(RANGE_C)];
        initial_sys = idgrey(@RCLadderGreyBox, initial_theta, 'c', [SCALING_R SCALING_C]);
        initial_sys.Structure.Parameters.Minimum = [repmat(RANGE_R(1), n, 1); repmat(RANGE_C(1), n, 1)];
        initial_sys.Structure.Parameters.Maximum = [repmat(RANGE_R(2), n, 1); repmat(RANGE_C(2), n, 1)];
        opt = greyestOptions('SearchMethod', 'auto', 'Focus', 'simulation');
        %opt.SearchOptions.Tolerance = 0.01;
        opt.SearchOptions.MaxIterations = 50;
        tic
        s(index).sys_idgrey = greyest(freq_id, initial_sys, opt);
        s(index).duration_idgrey = toc;
        s(index).sys_idgrey.Report.Termination
        s(index).theta_idgrey = getpvec(s(index).sys_idgrey).*[repmat(SCALING_R, n, 1); repmat(SCALING_C, n, 1)];
        fprintf('\nIDGREY: %i / %i correct', sum((abs(s(index).theta_idgrey-[s(index).R_true; s(index).C_true]))./[s(index).R_true; s(index).C_true]<0.01), 2*n)
        
        %% System Identification using N4SID and Structured Identification Algorithm
        tic
        s(index).sys_n4sid = n4sid(freq_id, n, 'FeedThrough', true);
        try
            s(index).sys_n4sid_est = RCLadderStructuredID(sys_est, 1);
            s(index).duration_n4sid_est = toc;
            [theta_R, theta_C] = RCLadder2Theta(s(index).sys_n4sid_est);
            s(index).theta_n4sid_est = [theta_R; theta_C];
        catch
            toc;
            s(index).sys_n4sid_est = [];
            s(index).duration_n4sid_est = nan;
            s(index).theta_n4sid_est = nan*ones(2*n, 1);
        end
    end
end
save(['IDMonteCarlo_Results_' datestr(now(), 30)], 's', '-v7.3')
fprintf('\n Finished\n.')

function distance = MatrixDistance(A, B)
    distance = norm(B - A)/norm(A);
end

function [A, B, C, D] = RCLadderGreyBox(theta, Ts, scaling)
    n = length(theta)/2;
    sys = RCLadderN(theta(1:n)*scaling(1), theta(n+1:end)*scaling(2), 'ascending', false);
    A = sys.A;
    B = sys.B;
    C = sys.C;
    D = sys.D;
    if Ts > 0
        A = eye(size(A)) + A * Ts; % Forward Euler discrete matrix
        B = Ts * B; % Forward Euler discrete matrix
    end
end