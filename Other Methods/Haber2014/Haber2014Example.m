%% Parameters from Section V
N = 40;
t_s = 1;
t_end = 1000;
A_true_local = [0.5728 0.1068; 0.1068 0.5728];
E_true_local = [0.1068 0; 0 0.1068];
B_true_local = [0.2136; 0.1068];
C_true_local = [1 0];
n = size(A_true_local, 1);
m = size(B_true_local, 2);
r = size(C_true_local, 1);

%% Build global system from (1)
A_true_global = kron(eye(N), A_true_local) + kron(diag(ones(N-1,1), 1), E_true_local) + kron(diag(ones(N-1,1), -1), E_true_local);
B_true_global = kron(eye(N), B_true_local);
C_true_global = kron(eye(N), C_true_local);
sys_true_global = ss(A_true_global, B_true_global, C_true_global, [], t_s);

%% Perform subspace identification of global system
t = 0:1/t_s:((t_end/t_s)-t_s);
u = randn(round(t_end/t_s), N);
y = lsim(sys_true_global, u, t);
data = iddata(y, u, t_s);
opt_n4sid = n4sidOptions('N4Horizon', [10, 15, 15], 'Focus', 'simulation');
sys_est_global = n4sid(data, N*n, 'Ts', t_s, opt_n4sid, 'DisturbanceModel', 'none')

%% Set up least-squares problem from (26)
[y_est, ~, x_est] = lsim(sys_est_global, u, t);
ii = 2;
k = 100;
%[x_est(k+1,((ii-1)*n+1):(ii*2))' [zeros(n-r, 1); y_est(k,1:r)']]
eq26_lhs = [x_est(k+1,3:4)' [zeros(n-r, 1); y_est(k,2)']];
eq26_rhs = [x_est(k,3:4)' zeros(2,1); x_est(k,1:2)' zeros(2,1); x_est(k,5:6)' zeros(2,1); u(k, 2) 0; zeros(2,1) x_est(k,3:4)'];
eq26_lhs/eq26_rhs