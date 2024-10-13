clear all 
close all
clc

% Inputs
n = 31;
nt = 1500;
x = linspace(0,3,n);
dx = x(2) - x(1);
gamma = 1.4;
c = 0.5;
throat = (n-1)/2;
a = 1 + 2.2*(x - 1.5).^2; % Area

[ sim_time_c, th_mach_c, th_press_c, th_temp_c, th_rho_c, mach_no_c, t_c, rho_c, p_c,m_dot_c ] = conservative( n, nt, x, dx, c, a, gamma, throat )

[ sim_time_nc, th_mach_nc, th_press_nc, th_temp_nc, th_rho_nc, mach_no_nc, t_nc, rho_nc, p_nc,m_dot_nc ] = non_conservative( n, nt, x, dx, c, a, gamma, throat )

fprintf('Simulation time for Conservative Form is %0.3g seconds', sim_time_c)
fprintf('\nSimulation time for Non-Conservative Form is %0.3g seconds', sim_time_nc)

% Plots
% Timewise Variation at the Nozzle Throat
figure(3)
subplot(411)
hold on
plot(th_mach_c, 'b')
plot(th_mach_nc, 'r')
legend('Conservative Form', 'Non-Conservative Form');
ylabel('Mach Number')
title('Timewise Variation at the Nozzle Throat')
grid minor

subplot(412)
hold on
plot(th_press_c, 'b')
plot(th_press_nc, 'r')
legend('Conservative Form', 'Non-Conservative Form');
ylabel('Pressure Ratio')
grid minor

subplot(413)
hold on
plot(th_rho_c, 'b')
plot(th_rho_nc, 'r')
legend('Conservative Form', 'Non-Conservative Form');
ylabel('Density Ratio')
grid minor

subplot(414)
hold on
plot(th_temp_c, 'b')
plot(th_temp_nc, 'r')
legend('Conservative Form', 'Non-Conservative Form');
xlabel('Number of Iterations')
ylabel('Temperature Ratio')
grid minor

% Qualitative Aspects of Quasi 1-D Nozzle Flow
figure(4)
subplot(411)
hold on
plot(x, mach_no_c, '-b+')
plot(x, mach_no_nc, 'r')
leg5 = legend('Conservative Form', 'Non-Conservative Form');
set(leg5, 'Location', 'northeastoutside')
ylabel('Mach Number')
title('Qualitative Aspects of Quasi 1-D Nozzle Flow')
grid minor

subplot(412)
hold on
plot(x, p_c, '-b+')
plot(x, p_nc, 'r')
leg6 = legend('Conservative Form', 'Non-Conservative Form');
set(leg6, 'Location', 'northeastoutside')
ylabel('Pressure Ratio')
grid minor

subplot(413)
hold on
plot(x, rho_c, '-b+')
plot(x, rho_nc, 'r')
leg7 = legend('Conservative Form', 'Non-Conservative Form');
set(leg7, 'Location', 'northeastoutside')
ylabel('Density Ratio')
grid minor

subplot(414)
hold on
plot(x, t_c, '-b+')
plot(x, t_nc, 'r')
leg8 = legend('Conservative Form', 'Non-Conservative Form');
set(leg8, 'Location', 'northeastoutside')
xlabel('Non-Dimensional Length of Nozzle')
ylabel('Temperature Ratio')
grid minor

figure(5)
hold on
plot(x, m_dot_c, 'b')
plot(x, m_dot_nc, 'r')
line([0 3], [0.579 0.579], 'color', 'g')
ylim([0.575 0.597])
leg8 = legend('Conservative Form', 'Non-Conservative Form', 'Exact Solution');
set(leg8, 'Location', 'northeastoutside')
title('Non-Dimensional Mass Flow Distribution');
xlabel('Non-Dimensional Length of Nozzle (x/l)')
ylabel('Mass Flow Rate Ratio')
grid minor

