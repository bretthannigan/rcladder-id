RANGE_R = [10 50]*1e3;
RANGE_C = [100 500]*1e-12;
IS_PRESCALE = true;
order = 14:14;
n_trials = 1;
s = struct('order', num2cell(repelem(order', n_trials)));
rng(1, 'twister') % Seed RNG
for i_n=1:length(order)
    for i_trial=1:n_trials
        %% Generate and Simulate Random True State-Space System.
        index = (i_n - 1)*n_trials + i_trial;
        n = order(i_n);
        R = RANGE_R(1) + rand(n, 1)*diff(RANGE_R);
        C = RANGE_C(1) + rand(n, 1)*diff(RANGE_C);
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
        [sys_est, T_est] = RCLadderStructuredID(s(index).sys_true, 1);
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
        
%         %% Run IDGREY Command
%         tic
%         sys_idgrey = idgrey(@RCLadderIDGrey, 1+rand(2*n, 1), 'c');
%         theta_idgrey = getpvec(sys_idgrey);
%         s(index).duration_idgrey = toc;
%         s(index).sys_idgrey = sys_idgrey;
%         s(index).R_idgrey = theta_idgrey(1:n);
%         s(index).C_idgrey = theta_idgrey(n+1:end);
    end
end
% %sum(sum(error(:,2:end), 2)>0.1)/size(error, 1)
% errors = reshape(sum(results(:,2:end), 2), [n_trials length(order)]);
% figure
% plot(order, 1-sum(errors>1e-4)/n_trials, '-*')
% xlabel('System order')
% ylabel('Fraction correctly reconstructed')
% title('Performance on Randomly Permuted Systems')
% ylim([0 1])

function distance = MatrixDistance(A, B)
    distance = norm(B - A)/norm(A);
end

function [A, B, C, D] = RCLadderIDGrey(theta, Ts)
    n = length(theta)/2;
    sys = RCLadderN(theta(1:n), theta(n+1:end), 'ascending', false);
    A = sys.A;
    B = sys.B;
    C = sys.C;
    D = sys.D;
end