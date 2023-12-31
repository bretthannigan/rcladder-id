order = 6;
RANGE_R = [1 10];
RANGE_C = [1 10];
SCALING_R = 1e4;
SCALING_C = 1e-12;
rng(1, 'twister') % Seed RNG
R = RANGE_R(1) + rand(order, 1)*diff(RANGE_R);
C = RANGE_C(1) + rand(order, 1)*diff(RANGE_C);

sys_true = RCLadderN(R*SCALING_R, C*SCALING_C, 'ascending', false);

freq_range = [1e4 1e7];
n_freq = 128;
w_test = 2*pi*logspace(log10(freq_range(1)), log10(freq_range(2)), n_freq);
[mag, phase, ~] = bode(sys_true, w_test);
mag = squeeze(mag);
phase = squeeze(phase);
response = mag.*exp(1i*deg2rad(phase));
freq_id = idfrd(response, w_test, 0);

%% Run IDGREY Command
initial_theta = [RANGE_R(1) + rand(order, 1)*diff(RANGE_R); RANGE_C(1) + rand(order, 1)*diff(RANGE_C)];
initial_sys = idgrey(@RCLadderGreyBox, initial_theta, 'c', [SCALING_R SCALING_C]);
initial_sys.Structure.Parameters.Minimum = [repmat(RANGE_R(1), order, 1); repmat(RANGE_C(1), order, 1)];
initial_sys.Structure.Parameters.Maximum = [repmat(RANGE_R(2), order, 1); repmat(RANGE_C(2), order, 1)];
opt = greyestOptions('SearchMethod', 'auto', 'Focus', 'simulation');
opt.SearchOptions.Tolerance = 1e-8;
%opt.SearchOptions.StepTolerance = eps;
%opt.SearchOptions.MaxIterations = 50;
tic
sys_idgrey = greyest(freq_id, initial_sys, opt)
toc
theta_idgrey = getpvec(sys_idgrey);

figure
bode(sys_true)
hold on
bode(sys_idgrey)

%% Run SSEST Command Plus Reconstruction Algorithm
tic
sys_est = n4sid(freq_id, order, 'FeedThrough', true)
[sys_est_structured, ~] = RCLadderStructuredID(sys_est, 1)
%sys_est = ssest(freq_id, idss(RCLadderN(initial_theta(1:n)*SCALING_R, initial_theta(n+1:end)*SCALING_C)))
toc
bode(sys_est)

function [A, B, C, D] = RCLadderGreyBox(theta, Ts, scaling)
    n = length(theta)/2;
    sys = RCLadderN(theta(1:n)*scaling(1), theta(n+1:end)*scaling(2), 'ascending', false);
    A = sys.A;
    B = sys.B;
    C = sys.C;
    D = sys.D;
    if Ts > 0 %sample the model with the sample time Ts
        A = eye(size(A)) + A * Ts; % Forward Euler discrete matrix
        B = Ts * B; % Forward Euler discrete matrix
    end
end