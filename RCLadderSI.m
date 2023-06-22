%% Create a true system with random R, C values.
n_section = 4;
rng(0, 'twister');
% Random non-negative resistance values from N(50k, 5k) distribution.
rvec = abs(10e3*randn(n_section, 1) + 100e3);
% Random non-negative capacitance values from N(50p, 10p) distribution.
cvec = abs(100e-12*randn(n_section, 1) + 10e-12);
sys_true = RCLadderN(rvec, cvec);
bode(sys_true)

%% Generate frequency response data.
t_s = 1/1e7;
f_start = 1e3; % Hz
f_end = 100e3; % Hz
n_f = 4;
f_test = logspace(log10(f_start), log10(f_end), n_f);
w_test = 2*pi*f_test;
sys_true_iddata = idfrd(sys_true, w_test);
h = squeeze(sys_true_iddata.ResponseData);

sysd_true = c2d(sys_true, t_s);
figure
bode(sys_true);
hold on
bode(sysd_true)
ax = get(gcf, 'Children');
ax = ax(2:end);
plot(ax(1), w_test, angle(h)*180/pi, '*')
plot(ax(2), w_test, 20*log10(abs(h)), '*')

sys_est = ssest(sys_true_iddata, n_section, 'Feedthrough', 1);

A_true = sys_true.A;
A_est = sys_est.A;

[L_true, U_true, P_true] = lu(inv(A_true));
LUIdentification2(A_true, rvec, cvec);

function LUIdentification1(A, r, c)
    [L, U, P] = lu(inv(A));
    n_A = length(A);
    X = L';
    X = [X(:,n_A)'; zeros(n_A-1, 1) triu(X(1:n_A-1, 2:n_A))];
    Y = X.*U;
    Y = Y + triu(Y, 1)';
    [U, S, V] = svd(-Y);
    
end