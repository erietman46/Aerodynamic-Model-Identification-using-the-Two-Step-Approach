%% A simple data file loader
% Coen de Visser, TU-Delft, 2026

clear; clc;

% Select datafile to load
[file, path] = uigetfile('*.mat');

if isequal(file,0)
    error('No file selected.');
end

dataname = fullfile(path, file);

try
    load(dataname);
catch
    error('Error loading data file');
end
%%
close all;
% variables in data file:
% t,alpha,Ax,Ay,Az,beta,da,de,dr,dta,dte,dtr,flaps,gamma,gear,Mach,p,phi,psi,q,r,Tc1,Tc2,theta,u_n,v_n,vtas,w_n
% string: contains aircraft parameters + variable explanation
figure(99);
hold on;
plot(t, p, 'b');
plot(t, q, 'r');
plot(t, r, 'k');
xlabel('time [s]');
legend('p [rad/s]', 'q [rad/s]', 'r [rad/s]');

figure(100);
hold on;
plot(t, phi, 'b');
plot(t, theta, 'r');
plot(t, psi, 'k');
xlabel('time [s]');
legend('\phi [rad]', '\theta [rad]', '\psi [rad]');

figure(101);
hold on;
plot(t, Ax, 'b');
plot(t, Ay, 'r');
plot(t, Az, 'k');
xlabel('time [s]');
legend('Ax [m/s^2]', 'Ay [m/s^2]', 'Az [m/s^2]');

figure(102);
hold on;
plot(t, alpha, 'b')
plot(t, beta, 'r');
xlabel('time [s]');
legend('alpha [rad]', 'beta [rad]');

figure(103);
hold on;
plot(t, vtas, 'b')
xlabel('time [s]');
legend('VTAS [m/s]');

figure(104);
hold on;
plot(t, Tc1, 'b')
plot(t, Tc2, 'r--')
xlabel('time [s]');
legend('Throttle (left)', 'Throttle (right)');

figure(105);
hold on;
plot(t, da, 'b')
plot(t, de, 'r');
plot(t, dr, 'k');
xlabel('time [s]');
legend('da [rad]', 'de [rad]', 'dr [rad]');
