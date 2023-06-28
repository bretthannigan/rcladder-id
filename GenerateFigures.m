%GENERATEFIGURES 
%   Script to generate the figures in Section XXX of [1]
%
%   See also: IDMONTECARLO
%       [1] Hannigan, B. C., Menon, C. (Draft). Fast, analytical method for 
%           structured identification of SISO RC-ladder-type systems.
%
%   $Author: BH$    $Date: 2023-06-28$  $Revision: 0$
%
%   ©2023 ETH Zurich, Brett Hannigan; D-HEST; Biomedical and Mobile Health Technology (BMHT) Lab; Carlo Menon

load('IDMonteCarlo_Results_20230628T130100.mat')

order = vertcat(s.order);
unique_order = unique(order);
A_est_dist = reshape(vertcat(s.A_est_dist), [], length(unique_order));
B_est_dist = reshape(vertcat(s.B_est_dist), [], length(unique_order));
C_est_dist = reshape(vertcat(s.C_est_dist), [], length(unique_order));
D_est_dist = reshape(vertcat(s.D_est_dist), [], length(unique_order));
max_est_dist = max(A_est_dist, B_est_dist);
max_est_dist = max(max_est_dist, C_est_dist);

A_hwang_dist = reshape(vertcat(s.A_hwang_dist), [], length(unique_order));
B_hwang_dist = reshape(vertcat(s.B_hwang_dist), [], length(unique_order));
C_hwang_dist = reshape(vertcat(s.C_hwang_dist), [], length(unique_order));
D_hwang_dist = reshape(vertcat(s.D_hwang_dist), [], length(unique_order));
max_hwang_dist = max(A_hwang_dist, B_hwang_dist);
max_hwang_dist = max(max_hwang_dist, C_hwang_dist);

figure
plot(unique_order, sum(max_est_dist<0.01)/size(A_est_dist, 1), 'o-', 'DisplayName', 'This work')
hold on
plot(unique_order, sum(max_hwang_dist<0.01)/size(A_hwang_dist, 1), 'o-', 'DisplayName', 'Hwang1984')
xlabel('System order')
ylabel('Fraction correctly reconstucted')
legend('Location', 'NorthEast')

n_R_est_correct = zeros(length(unique_order), 1);
n_C_est_correct = zeros(length(unique_order), 1);
n_R_idgrey_correct = zeros(length(unique_order), 1);
n_C_idgrey_correct = zeros(length(unique_order), 1);

for i_order=1:length(unique_order)
    R_true = horzcat(s(order==unique_order(i_order)).R_true)';
    C_true = horzcat(s(order==unique_order(i_order)).C_true);
    theta_est = horzcat(s(order==unique_order(i_order)).theta_n4sid_est)';
    R_est = theta_est(:,1:unique_order(i_order));
    C_est = theta_est(:,unique_order(i_order)+1:end);
    n_R_est_correct(i_order) = mean(sum((abs(R_est-R_true)./R_true)<=0.01, 2));
    n_C_est_correct(i_order) = mean(sum((abs(C_est-C_true)./C_true)<=0.01 ,2));
    theta_idgrey = horzcat(s(order==unique_order(i_order)).theta_idgrey)';
    R_idgrey = theta_idgrey(:,1:unique_order(i_order));
    C_idgrey = theta_idgrey(:,unique_order(i_order)+1:end);
    n_R_idgrey_correct(i_order) = mean(sum((abs(R_idgrey-R_true)./R_true)<=0.01, 2));
    n_C_idgrey_correct(i_order) = mean(sum((abs(C_idgrey-C_true)./C_true)<=0.01 ,2));
end

figure
plot(unique_order, n_R_est_correct, '-o', 'DisplayName', 'R, this work')
hold on
plot(unique_order, n_C_est_correct, '-o', 'DisplayName', 'C, this work')
plot(unique_order, n_R_idgrey_correct, '-o', 'DisplayName', 'R, idgrey')
plot(unique_order, n_C_idgrey_correct, '-o', 'DisplayName', 'C, idgrey')
plot([unique_order(1) unique_order(end)], [unique_order(1) unique_order(end)], '--k', 'DisplayName', 'Perfect')
xlabel('System order')
ylabel('Number correct parameters')
legend('Location', 'NorthWest')