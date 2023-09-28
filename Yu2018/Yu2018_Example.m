%Test with an example from [Section III-A, Prot2020]
fn_A = @(theta) [-theta(1) theta(3) 0; theta(1) -theta(2)-theta(3) theta(4); 0 theta(2) -theta(4)];
fn_B = @(theta) [0; 0; 1];
fn_C = @(theta) [0 0 1];

theta = [-0.394 -0.893 0.325 0.383];

sys = ss(fn_A(theta), fn_B(theta), fn_C(theta), [], 1);
T = randn(3);
sys_t = ss2ss(sys, T);

[sys_est, theta_est] = Yu2018(sys_t, fn_A, fn_B, fn_C, length(theta));

disp(norm(theta_est - theta))

% A_star = [-0.3977 -0.1384 -0.0944 -0.2233 -0.1641;...
%           -0.1596 -0.3357 -0.1360 -0.1127 -0.2045;...
%           -0.0087 -0.0776 -0.1824 -0.0326 -0.2518;...
%           -0.2079 -0.1306  0.0821 -0.3278  0.0146;...
%           -0.1899 -0.2372 -0.1908  0.0601 -0.6153];
% B_star = [0.9509; -1.2930; 0.4921; -0.9242; 1.0427];
% C_star = [0.1037 -0.3812 0.5442 -0.9107 2.4788];
% n = size(A_star, 1);
% n_omit = 1;
% i_omit = randi(n, n_omit, 2);
% 
% fn_A = @(theta) [-0.3977 -0.1384 -0.0944 -0.2233 -0.1641;...
%                  -0.1596 -0.3357 -0.1360 -0.1127 -0.2045;...
%                  -0.0087 -0.0776 -0.1824 theta(1) -0.2518;...
%                  -0.2079 -0.1306  0.0821 -0.3278  0.0146;...
%                  -0.1899 -0.2372 -0.1908  0.0601 -0.6153];
% fn_B = @(theta) B_star;
% fn_C = @(theta) C_star;
% 
% sys = ss(A_star, B_star, C_star, [], 1);
% theta_est = Yu2018(sys, fn_A, fn_B, fn_C, n_omit)