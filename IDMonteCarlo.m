RANGE_R = [1 10];
RANGE_C = [1 10];
SCALING_R = 1e4;
SCALING_C = 1e-12;
IS_PRESCALE = true;
order = 4:10;
n_trials = 100;
s = struct('order', num2cell(repelem(order', n_trials)));
rng(1, 'twister') % Seed RNG
for i_n=1:length(order)
    for i_trial=1:n_trials
        %% Generate and Simulate Random True State-Space System.
        index = (i_n - 1)*n_trials + i_trial;
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
        tic
        [sys_est, T_est] = RCLadderStructuredID(s(index).sys_t, 1);
        s(index).duration_est = toc;
        s(index).sys_est = sys_est;
        s(index).A_est_dist = MatrixDistance(s(index).sys_true.A, s(index).sys_est.A);
        s(index).B_est_dist = MatrixDistance(s(index).sys_true.B, s(index).sys_est.B);
        s(index).C_est_dist = MatrixDistance(s(index).sys_true.C, s(index).sys_est.C);
        s(index).D_est_dist = MatrixDistance(s(index).sys_true.D, s(index).sys_est.D);
        
        %% Run Routh Array Procedure
        tic
        [sys_hwang, ~] = RCLadderHwang(s(index).sys_t);
        s(index).duration_hwang = toc;
        s(index).sys_hwang = sys_hwang;
        s(index).A_hwang_dist = MatrixDistance(s(index).sys_true.A, s(index).sys_hwang.A);
        s(index).B_hwang_dist = MatrixDistance(s(index).sys_true.B, s(index).sys_hwang.B);
        s(index).C_hwang_dist = MatrixDistance(s(index).sys_true.C, s(index).sys_hwang.C);
        s(index).D_hwang_dist = MatrixDistance(s(index).sys_true.D, s(index).sys_hwang.D);

        
    end
end

function distance = MatrixDistance(A, B)
    distance = norm(B - A)/norm(A);
end