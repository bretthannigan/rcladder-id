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

load('IDMonteCarlo_Results_20230703T202116.mat')

COMPARISON_THRESHOLD = 0.01;

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

set(0,'DefaultFigureColor','remove')
fig1 = figure();
ax(1) = subplot(2, 2, 1);
set(ax(1), 'pos', [0.07 0.57 0.4 0.4]);
plot(unique_order, sum(max_est_dist<0.01)/size(A_est_dist, 1), 'o-k', 'LineWidth',  1, 'MarkerFaceColor', 'white', 'DisplayName', 'This work')
hold on
plot(unique_order, sum(max_hwang_dist<0.01)/size(A_hwang_dist, 1), 'o-', 'Color', [0.4 0.4 0.4], 'LineWidth',  1, 'MarkerFaceColor', 'white', 'DisplayName', 'Ref. [9]')
xlabel('System order')
ylabel('Fraction correctly reconstucted')
legend('Location', 'SouthWest')
xlim([unique_order(1) unique_order(end)])

n_R_est_correct = zeros(length(unique_order), 1);
n_C_est_correct = zeros(length(unique_order), 1);
n_R_idgrey_correct = zeros(length(unique_order), 1);
n_C_idgrey_correct = zeros(length(unique_order), 1);

for i_order=1:length(unique_order)
    R_true = horzcat(s(order==unique_order(i_order)).R_true)';
    C_true = horzcat(s(order==unique_order(i_order)).C_true)';
    theta_est = horzcat(s(order==unique_order(i_order)).theta_n4sid_est)';
    R_est = theta_est(:,1:unique_order(i_order));
    C_est = theta_est(:,unique_order(i_order)+1:end);
    n_R_est_correct(i_order) = mean(sum((abs(R_est-R_true)./R_true)<=COMPARISON_THRESHOLD, 2));
    n_C_est_correct(i_order) = mean(sum((abs(C_est-C_true)./C_true)<=COMPARISON_THRESHOLD ,2));
    theta_idgrey = horzcat(s(order==unique_order(i_order)).theta_idgrey)';
    R_idgrey = theta_idgrey(:,1:unique_order(i_order));
    C_idgrey = theta_idgrey(:,unique_order(i_order)+1:end);
    n_R_idgrey_correct(i_order) = mean(sum((abs(R_idgrey-R_true)./R_true)<=COMPARISON_THRESHOLD, 2));
    n_C_idgrey_correct(i_order) = mean(sum((abs(C_idgrey-C_true)./C_true)<=COMPARISON_THRESHOLD ,2));
end

ax(2) = subplot(2, 2, 2);
set(ax(2), 'pos', [0.55 0.57 0.4 0.4]);
plot([0  unique_order(end)], [0 unique_order(end)], '--', 'Color', [0.6 0.6 0.6], 'DisplayName', 'Perfect')
hold on
plot(unique_order, n_R_est_correct, '-squarek', 'MarkerSize', 8, 'LineWidth',  1, 'MarkerFaceColor', 'white', 'DisplayName', 'R, this work')
plot(unique_order, n_C_est_correct, '-^k', 'LineWidth',  1, 'MarkerFaceColor', 'white', 'DisplayName', 'C, this work')
plot(unique_order, n_R_idgrey_correct, '-square', 'Color', [0.4 0.4 0.4], 'MarkerSize', 8, 'LineWidth',  1, 'MarkerFaceColor', 'white', 'DisplayName', 'R, idgrey')
plot(unique_order, n_C_idgrey_correct, '-^', 'Color', [0.4 0.4 0.4], 'LineWidth',  1, 'MarkerFaceColor', 'white', 'DisplayName', 'C, idgrey')
xlabel('System order')
ylabel('Number of correct parameters')
legend('Location', 'NorthWest')
xlim([unique_order(1) unique_order(end)])

duration_est = nanmedian(reshape(vertcat(s.duration_est), [], length(unique_order)));
duration_hwang = nanmedian(reshape(vertcat(s.duration_hwang), [], length(unique_order)));
duration_n4sid_est = nanmedian(reshape(vertcat(s.duration_n4sid_est), [], length(unique_order)));
duration_idgrey = nanmedian(reshape(vertcat(s.duration_idgrey), [], length(unique_order)));

ax(3) = subplot(2, 2, 3);
set(ax(3), 'pos', [0.07 0.08 0.4 0.4]);
plot(unique_order, duration_est, '-ok',  'LineWidth', 1, 'MarkerFaceColor', 'white', 'DisplayName', 'this work')
hold on
plot(unique_order, duration_hwang, '-o', 'Color', [0.4 0.4 0.4], 'LineWidth',  1, 'MarkerFaceColor', 'white', 'DisplayName', 'Ref. [9]')
xlabel('System order')
ylabel('Computation time (s)')
legend('Location', 'NorthWest')
xlim([unique_order(1) unique_order(end)])

ax(4) = subplot(2, 2, 4);
set(ax(4), 'pos', [0.55 0.08 0.4 0.4]);
plot(unique_order, duration_n4sid_est, '-ok', 'LineWidth',  1, 'MarkerFaceColor', 'white', 'DisplayName', 'this work')
hold on
plot(unique_order, duration_idgrey, '-o', 'Color', [0.4 0.4 0.4], 'LineWidth',  1, 'MarkerFaceColor', 'white', 'DisplayName', 'idgrey')
xlabel('System order')
ylabel('Computation time (s)')
legend('Location', 'NorthWest')
xlim([unique_order(1) unique_order(end)])

fig1.Units               = 'centimeters';
fig1.Position(3)         = 18;
fig1.Position(4)         = 12;
set(fig1.Children, 'FontName', 'Times', 'FontSize', 9);
fig1.PaperPositionMode = 'auto';
print('img/Comparison-Figure.png', '-dpng', '-r600')
print('img/Comparison-Figure.eps', '-depsc')