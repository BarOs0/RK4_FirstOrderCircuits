% ==================================================================================
% AGH University of Science and Technology
% Authors: Michał Dobrzycki, Bartłomiej Osiak
% Faculty: WIEiT
% Field: EiT, II year, winter semester
% Project: Frequency Analysis
% Date: 27th January 2025
% ==================================================================================
% Description:
% This script performs frequency analysis using the fourth-order Runge-Kutta
% method for numerical integration. The analysis is applied to RL and RC circuits 
% driven by a sine-wave input signal.
% ===================================================================================
%
%            CR        RC
%          o   o     o   o
%          |   |     |   |
%       o--o[R]o-----o|C|o--+
%                           |
%      (~)                  |
%                           |
%       o-------------------+
%
%            RL        LR
%          o   o     o   o
%          |   |     |   |
%       o--o(L)o-----o[R]o--+
%                           |
%      (~)                  |
%                           |
%       o-------------------+
%
% =====================================================================================
%% Initialization

clear; clc; close all

% Time Configurations
tstart = 0; % start time [s]  
h = 5e-6; % time step [s]
tend = 2; % end time [s]
t = tstart:h:tend; % time vector

% Frequency vector
freq = logspace(0.01, 4, 50);

% Input wave configurations
A = 50*1e-3; % Amplitude of the input signal [V]
F = 300; % Frequency of the input signal [Hz]
Vin = Fnk(A, F, t); % Input voltage (sinusoidal)

% Circuit Configurations
R = 1e3; % [Ohms]
C = 1e-6; % [Farads]
L = 1; % [Henrys]

% Initial conditions
y0 = 0;  

% Analytical Bode Characteristics 
H_RC = Transmitance(R,L,C,freq,'RC');
H_CR = Transmitance(R,L,C,freq,'CR');
H_LR = Transmitance(R,L,C,freq,'LR');
H_RL = Transmitance(R,L,C,freq,'RL');

%% Computing values for RC, CR, LR, RL
vpp_RC = []; 
vpp_CR = [];
phase_RC = [];
phase_CR = [];
y1 = RK4(t, y0, R, C, L, h, A, F, 'c'); 
y2 = Vin-y1; 

for k = freq
    Vin_k = Fnk(A, k, t);          
    y_c = RK4(t, y0, R, C, L, h, A, k, 'c');
    phase_RC(end+1) = ComputePhaseShift(Vin_k, y_c, k, t);
    phase_CR(end+1) = ComputePhaseShift(Vin_k, Vin_k-y_c, k, t);
    vpp_RC(end+1) = ComputeVpp(y_c);
    vpp_CR(end+1) = ComputeVpp(Vin_k-y_c);
end
gain1 = vpp_RC./(2*A); 
gain2 = vpp_CR./(2*A);

vpp_LR = []; 
vpp_RL = [];
phase_LR = [];
phase_RL = [];
y3 = RK4(t, y0, R, C, L, h, A, F, 'l').*R; 
y4 = Vin-y3;

for k = freq
    Vin_k = Fnk(A, k, t);          
    y_l = RK4(t, y0, R, C, L, h, A, k, 'l').*R;
    phase_LR(end+1) = ComputePhaseShift(Vin_k, y_l, k, t);
    phase_RL(end+1) = ComputePhaseShift(Vin_k, Vin_k-y_l, k, t);
    vpp_LR(end+1) = ComputeVpp(y_l);       
    vpp_RL(end+1) = ComputeVpp(Vin_k-y_l);
end
gain3 = vpp_LR./(2*A); 
gain4 = vpp_RL./(2*A);

%% Output response plots
figure(1);
sgtitle('Comparing output to the input');

subplot(2,2,1); 
    plot(t, y1, 'Color', '#ffb266', 'LineWidth', 1.5); hold on;
    plot(t, Vin, 'b', 'LineWidth', 1.5); grid on;
    xlim([1,1.01]); % Setting precision for time dimension
    title('RC');
    xlabel('t [s]');
    ylabel('U [V]');
    legend('V_{out}');
    
subplot(2,2,2);
    plot(t, y2, 'c', 'LineWidth', 1.5); hold on;
    plot(t, Vin, 'b', 'LineWidth', 1.5); grid on;
    xlim([1,1.01]);
    title('CR');
    xlabel('t [s]');
    ylabel('U [V]');
    legend('V_{out}');
    
subplot(2,2,3);
    plot(t, y3, 'g', 'LineWidth', 1.5); hold on;
    plot(t, Vin, 'b', 'LineWidth', 1.5); grid on;
    xlim([1,1.01]);
    title('LR');
    xlabel('t [s]');
    ylabel('U [V]');
    legend('V_{out}');
    
subplot(2,2,4);
    plot(t, y4, 'm', 'LineWidth', 1.5); hold on;
    plot(t, Vin, 'b', 'LineWidth', 1.5); grid on;
    xlim([1,1.01]);
    title('RL');
    xlabel('t [s]');
    ylabel('U [V]');
    legend('V_{out}');

%% Frequency plots
f3db_RC_CR = (1/(R*C))/(2*pi); % Calculating analytical poles
f3db_RL_LR = (R/L)/(2*pi);

figure(2)
sgtitle('RC');
subplot(2,1,1)
semilogx(freq, 20*log10(gain1), 'Color', '#ffb266', 'LineWidth', 2); hold on;
semilogx(freq, 20*log10(abs(H_RC)), 'k.', 'LineWidth', 2); hold on;
semilogx(f3db_RC_CR, -3, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k'); hold on; 
xline(f3db_RC_CR, '--k', 'LineWidth', 1);
grid on;
legend('Amplitude', 'Analytical Amplitude', 'f_{3dB}', 'Location', 'SouthWest');
xlabel('Frequency [Hz]');
ylabel('Gain [dB]');
ylim([-45,5]);
subplot(2,1,2)
semilogx(freq, phase_RC, 'Color', '#ffb266', 'LineWidth', 2); hold on;
semilogx(freq, angle(H_RC)*(180/pi), 'k.', 'LineWidth', 2); grid on;
legend('Phase', 'Analytical Phase', 'Location', 'SouthWest');
xlabel('Frequency [Hz]');
ylabel('Phase [degrees]');

figure(3)
sgtitle('CR');
subplot(2,1,1)
semilogx(freq, 20*log10(gain2), 'c', 'LineWidth', 2); hold on;
semilogx(freq, 20*log10(abs(H_CR)), 'k.', 'LineWidth', 2); hold on;
semilogx(f3db_RC_CR, -3, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k'); hold on; 
xline(f3db_RC_CR, '--k', 'LineWidth', 1);
grid on;
legend('Amplitude', 'Analytical Amplitude', 'f_{3dB}', 'Location', 'SouthEast');
xlabel('Frequency [Hz]');
ylabel('Gain [dB]');
ylim([-45,5]);
subplot(2,1,2)
semilogx(freq, phase_CR, 'c', 'LineWidth', 2); hold on;
semilogx(freq, angle(H_CR)*(180/pi), 'k.', 'LineWidth', 2); grid on;
legend('Phase', 'Analytical Phase', 'Location', 'SouthWest');
xlabel('Frequency [Hz]');
ylabel('Phase [degrees]');

figure(4)
sgtitle('LR');
subplot(2,1,1)
semilogx(freq, 20*log10(gain3), 'g', 'LineWidth', 2); hold on;
semilogx(freq, 20*log10(abs(H_LR)), 'k.', 'LineWidth', 2); hold on;
semilogx(f3db_RL_LR, -3, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k'); hold on; 
xline(f3db_RL_LR, '--k', 'LineWidth', 1);
grid on;
legend('Amplitude', 'Analytical Amplitude', 'f_{3dB}', 'Location', 'SouthWest');
xlabel('Frequency [Hz]');
ylabel('Gain [dB]');
ylim([-45,5]);
subplot(2,1,2)
semilogx(freq, phase_LR, 'g', 'LineWidth', 2); hold on;
semilogx(freq, angle(H_LR)*(180/pi), 'k.', 'LineWidth', 2); grid on;
legend('Phase', 'Analytical Phase', 'Location', 'SouthWest');
xlabel('Frequency [Hz]');
ylabel('Phase [degrees]');

figure(5)
sgtitle('RL');
subplot(2,1,1)
semilogx(freq, 20*log10(gain4), 'm', 'LineWidth', 2); hold on;
semilogx(freq, 20*log10(abs(H_RL)), 'k.', 'LineWidth', 2); hold on;
semilogx(f3db_RL_LR, -3, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k'); hold on; 
xline(f3db_RL_LR, '--k', 'LineWidth', 1);
grid on;
legend('Amplitude', 'Analytical Amplitude', 'f_{3dB}', 'Location', 'SouthEast');
xlabel('Frequency [Hz]');
ylabel('Gain [dB]');
ylim([-45,5]);
subplot(2,1,2)
semilogx(freq, phase_RL, 'm', 'LineWidth', 2); hold on;
semilogx(freq, angle(H_RL)*(180/pi), 'k.', 'LineWidth', 2); grid on;
legend('Phase', 'Analytical Phase', 'Location', 'SouthWest');
xlabel('Frequency [Hz]');
ylabel('Phase [degrees]');
