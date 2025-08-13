% Script that generates Figures 3, 4, and 5 in "Improving frequency
% stability using slowly modulated adaptive feedback," by H. Dankowicz,
% S.W. Shaw, and O. Shoshani.

% The code contains five simulations (numbered from 1 to 5) of three
% different cases: 
%   1) Standard Duffing Oscillator (stops after 100 time units). 
%   2) Oscillator with zero amplitude-to-frequency noise conversion (stops 
%       after 100 time units). 
%   3) Oscillator with theoretically zero phase diffusion (stops after 
%       100 time units). 
%   4) Oscillator with theoretically zero phase diffusion (stops after 
%       50 time units).
%   5) Oscillator with theoretically zero phase diffusion (stops after 
%       200 time units).

% Figure 5 of the manuscript can also be generated with this code by
% defining D_phi = sigma^2 and using different values of That and k.

%% 1) The Standard Duffing Oscillator

% System parameters:
Gamma1   = 1;
gamma    = 1;
Delta    = pi/2;
omega1   = 1;
S        = 3;
k        = 0;

% Initial state
That     = 0;
thetahat = 0;
a1hat    = (S*sin(Delta)-That*sin(thetahat))/2/omega1/Gamma1;
phi1hatp = (3*gamma*a1hat^3-4*S*cos(Delta)+4*That*cos(thetahat))/8/omega1/a1hat;

% Noise intensities
Ixi  = 1e-4;
Ieta = 1e-4;

M  = 1000; %Size of ensemble
dt = 10^(-4);
N  = 1000000; %time length
t  = (0:dt:N*dt)';

dist_duff = zeros(M,1);
phi       = zeros(M,N+1);
for i=1:M
    a1    = a1hat;
    T     = That;
    theta = thetahat;
    phi1  = phi1hatp;
    Xem = [a1 T theta phi1; zeros(N,4)];

    for j=1:N
        dW = sqrt(dt)*randn(1,2);
        Ia11   = (omega1^2)*(a1^2)*Ieta/4 + Ixi/(omega1^2);
        Iphi11 = 3*(omega1^2)*Ieta/4 + Ixi/((a1^2)*(omega1^2));
        Fnoise = [sqrt(Ia11)*dW(1); sqrt(Iphi11)*dW(2)];
        f = [-Gamma1*a1+(S*sin(Delta)-T*sin(theta))/2/omega1; ...
            a1-a1hat+That-T; ...
            a1-a1hat+2*T-2*That; ...
            (3*gamma*a1^3-4*S*cos(Delta)+4*T*cos(theta))/8/omega1/a1];
        a1    = a1    + dt*f(1) + Fnoise(1);
        T     = 0;
        theta = theta + dt*f(3) + k*Fnoise(2);
        phi1  = phi1  + dt*f(4) + Fnoise(2);
        Xem(j+1,:) = [a1 T theta phi1];
    end
    phi(i,:) = Xem(:,4)';
    dist_duff(i) = Xem(end,4);
end

DT_duff = Iphi11+Ia11*(3*gamma*a1hat/(4*omega1*Gamma1))^2;
dp_duff = var(phi);
y1 = -0.1:0.0001:0.1;
mu1 = 0;
sigma1 = sqrt(var(dist_duff/sqrt(100)));
f1 = exp(-(y1-mu1).^2./(2*sigma1^2))./(sigma1*sqrt(2*pi));

figure(1)
plot(t, dp_duff, 'b', 'LineWidth', 2)
hold on
plot(t, DT_duff*t, '--b', 'LineWidth', 2);
set(gca, 'FontSize', 22)
xlabel('$t$', 'Interpreter', 'Latex', 'Fontsize', 24)
ylabel('$\langle\phi^2(s)\rangle-\langle\phi(s)\rangle^2$', ...
    'Interpreter', 'Latex', 'Fontsize', 24)

figure(2)
histogram(dist_duff/sqrt(100)-mean(dist_duff/sqrt(100)), ...
    'LineWidth', 2, 'EdgeColor', [1 1 1], 'FaceColor', [0 0 1], ...
    'Normalization', 'pdf', 'NumBins', 10)
hold on
plot(y1, f1, 'b--', 'LineWidth', 2)
set(gca, 'FontSize', 22)
ylabel('$w(\delta\phi(s)/\sqrt{s})$', 'Interpreter', 'Latex', ...
    'Fontsize', 24)
xlabel('$\delta\phi(s)/\sqrt{s}$', 'Interpreter', 'Latex', ...
    'Fontsize', 24)

%% 2) Oscillator with zero A-f noise conversion

clear all

Gamma1   = 1;
gamma    = 1;
Delta    = pi/2;
omega1   = 1;
S        = 3;
k        = 0;

That     = 1;
thetahat = 0;
a1hat    = (S*sin(Delta)-That*sin(thetahat))/2/omega1/Gamma1;
phi1hatp = (3*gamma*a1hat^3-4*S*cos(Delta)+4*That*cos(thetahat))/8/omega1/a1hat;

Ixi=1e-4;
Ieta=1e-4;

M  = 1000;
dt = 10^(-4);
N  = 1000000;
t  = (0:dt:N*dt)';

dist_nAF = zeros(M,1);
phi      = zeros(M,N+1);
for i=1:M
    a1    = a1hat;
    T     = That;
    theta = thetahat;
    phi1  = phi1hatp;
    Xem = [a1 T theta phi1; zeros(N,4)];

    for j=1:N
        dW = sqrt(dt)*randn(1,2);
        Ia11   = (omega1^2)*(a1^2)*Ieta/4 + Ixi/(omega1^2);
        Iphi11 = 3*(omega1^2)*Ieta/4 + Ixi/((a1^2)*(omega1^2));
        Fnoise = [sqrt(Ia11)*dW(1); sqrt(Iphi11)*dW(2)];
        f = [-Gamma1*a1+(S*sin(Delta)-T*sin(theta))/2/omega1; ...
            a1-a1hat+That-T; ...
            a1-a1hat+2*T-2*That; ...
            (3*gamma*a1^3-4*S*cos(Delta)+4*T*cos(theta))/8/omega1/a1];
        a1    = a1    + dt*f(1) + Fnoise(1);
        T     = T     + dt*f(2);
        theta = theta + dt*f(3) + k*Fnoise(2);
        phi1  = phi1  + dt*f(4) + Fnoise(2);
        Xem(j+1,:) = [a1 T theta phi1];
    end
    phi(i,:) = Xem(:,4)';
    dist_nAf(i) = Xem(end,4);
end

DT_nAf = Iphi11;
dp_nAf = var(phi);
y2 = -0.1:0.0001:0.1;
mu2 = 0;
sigma2 = sqrt(var(dist_nAf/sqrt(100)));
f2 = exp(-(y2-mu2).^2./(2*sigma2^2))./(sigma2*sqrt(2*pi));

figure(1)
plot(t, dp_nAf, 'k', 'LineWidth', 2)
plot(t, DT_nAf*t, '--k', 'LineWidth', 2);

figure(3)
histogram(dist_nAf/sqrt(100)-mean(dist_nAf/sqrt(100)), ...
    'LineWidth', 2, 'EdgeColor', [1 1 1], 'FaceColor', [0 0 0], ...
    'Normalization', 'pdf', 'NumBins', 10)
hold on
plot(y2, f2, 'k--', 'LineWidth', 2)
set(gca, 'FontSize', 22)
ylabel('$w(\delta\phi(s)/\sqrt{s})$', 'Interpreter', 'Latex', ...
    'Fontsize', 24)
xlabel('$\delta\phi(s)/\sqrt{s}$', 'Interpreter', 'Latex', ...
    'Fontsize', 24)

%% 3) The zero diffusion Oscillator

clear all

Gamma1   = 1;
gamma    = 1;
Delta    = pi/2;
omega1   = 1;
S        = 3;
k        = 216/89;

That     = 1;
thetahat = 0;
a1hat    = (S*sin(Delta)-That*sin(thetahat))/2/omega1/Gamma1;
phi1hatp = (3*gamma*a1hat^3-4*S*cos(Delta)+4*That*cos(thetahat))/8/omega1/a1hat;

Ixi  = 1e-4;
Ieta = 1e-4;

M  = 1000;
dt = 10^(-4);
N  = 1000000;
t  = (0:dt:N*dt)';

dist_opt = zeros(M,1);
phi      = zeros(M,N+1);
for i=1:M
    a1    = a1hat;
    T     = That;
    theta = thetahat;
    phi1  = phi1hatp;
    Xem = [a1 T theta phi1; zeros(N,4)];

    for j=1:N
        dW = sqrt(dt)*randn(1,2);
        Ia11   =(omega1^2)*(a1^2)*Ieta/4 + Ixi/(omega1^2);
        Iphi11 = 3*(omega1^2)*Ieta/4 + Ixi/((a1^2)*(omega1^2));
        Fnoise = [sqrt(Ia11)*dW(1); sqrt(Iphi11)*dW(2)];
        f = [-Gamma1*a1+(S*sin(Delta)-T*sin(theta))/2/omega1; ...
            a1-a1hat+That-T; ...
            a1-a1hat+2*T-2*That; ...
            (3*gamma*a1^3-4*S*cos(Delta)+4*T*cos(theta))/8/omega1/a1];
        a1    = a1    + dt*f(1) + Fnoise(1);
        T     = T     + dt*f(2);
        theta = theta + dt*f(3) + k*Fnoise(2);
        phi1  = phi1  + dt*f(4) + Fnoise(2);
        Xem(j+1,:) = [a1 T theta phi1];
    end
    phi(i,:) = Xem(:,4)';
    dist_opt(i) = Xem(end,4);
end

DT_opt = 0;
dp_opt = var(phi);
y3 = -0.1:0.0001:0.1;
mu3 = 0;
sigma3 = sqrt(var(dist_opt/sqrt(100)));
f3 = exp(-(y3-mu3).^2./(2*sigma3^2))./(sigma3*sqrt(2*pi));

figure(1)
plot(t, dp_opt, 'r', 'LineWidth', 2)
plot(t, DT_opt*t, '--r', 'LineWidth', 2);

figure(4)
histogram(dist_opt/sqrt(100)-mean(dist_opt/sqrt(100)), ...
    'LineWidth', 0.5, 'EdgeColor', [1 1 1], 'FaceColor', [1 0 0], ...
    'Normalization', 'pdf', 'NumBins', 10)
hold on
plot(y3, f3, 'r--', 'LineWidth', 2)
set(gca, 'FontSize', 22)
ylabel('$w(\delta\phi(s)/\sqrt{s})$', 'Interpreter', 'Latex', ...
    'Fontsize', 24)
xlabel('$\delta\phi(s)/\sqrt{s}$', 'Interpreter', 'Latex', ...
    'Fontsize', 24)

%% 4) The zero diffusion Oscillator short

clear all

Gamma1   = 1;
gamma    = 1;
Delta    = pi/2;
omega1   = 1;
S        = 3;
k        = 216/89;

That     = 1;
thetahat = 0;
a1hat    = (S*sin(Delta)-That*sin(thetahat))/2/omega1/Gamma1;
phi1hatp = (3*gamma*a1hat^3-4*S*cos(Delta)+4*That*cos(thetahat))/8/omega1/a1hat;

Ixi  = 1e-4;
Ieta = 1e-4;

M  = 1000;
dt = 10^(-4);
N  = 500000;
t  = (0:dt:N*dt)';

dist_opt_s = zeros(M,1);
phi        = zeros(M,N+1);
for i=1:M
    a1    = a1hat;
    T     = That;
    theta = thetahat;
    phi1  = phi1hatp;
    Xem = [a1 T theta phi1; zeros(N,4)];

    for j=1:N
        dW = sqrt(dt)*randn(1,2);
        Ia11   = (omega1^2)*(a1^2)*Ieta/4 + Ixi/(omega1^2);
        Iphi11 = 3*(omega1^2)*Ieta/4 + Ixi/((a1^2)*(omega1^2));
        Fnoise = [sqrt(Ia11)*dW(1); sqrt(Iphi11)*dW(2)];
        f = [-Gamma1*a1+(S*sin(Delta)-T*sin(theta))/2/omega1; ...
            a1-a1hat+That-T; ...
            a1-a1hat+2*T-2*That; ...
            (3*gamma*a1^3-4*S*cos(Delta)+4*T*cos(theta))/8/omega1/a1];
        a1    = a1    + dt*f(1) + Fnoise(1);
        T     = T     + dt*f(2);
        theta = theta + dt*f(3) + k*Fnoise(2);
        phi1  = phi1  + dt*f(4) + Fnoise(2);
        Xem(j+1,:) = [a1 T theta phi1];
    end
    phi(i,:) = Xem(:,4)';
    dist_opt_s(i) = Xem(end,4);
end

y4 = -0.01:0.0001:0.01;
mu4 = 0;
sigma4 = sqrt(var(dist_opt_s/sqrt(50)));
f4 = exp(-(y4-mu4).^2./(2*sigma4^2))./(sigma4*sqrt(2*pi));

figure(5)
histogram(dist_opt_s/sqrt(100)-mean(dist_opt_s/sqrt(100)), ...
    'LineWidth', 2, 'EdgeColor', [1 1 1], 'FaceColor', [1 0 0], ...
    'Normalization', 'pdf', 'NumBins', 10)
hold on
plot(y4, f4, 'r--', 'LineWidth', 2)
xlim([-0.01,0.01])
set(gca, 'FontSize', 22)
ylabel('$w(\delta\phi(s)/\sqrt{s})$', 'Interpreter', 'Latex', ...
    'Fontsize', 24)
xlabel('$\delta\phi(s)/\sqrt{s}$', 'Interpreter', 'Latex', ...
    'Fontsize', 24)

%% 5) The zero diffusion Oscillator long

clear all

Gamma1   = 1;
gamma    = 1;
Delta    = pi/2;
omega1   = 1;
S        = 3;
k        = 216/89;

That     = 1;
thetahat = 0;
a1hat    = (S*sin(Delta)-That*sin(thetahat))/2/omega1/Gamma1;
phi1hatp = (3*gamma*a1hat^3-4*S*cos(Delta)+4*That*cos(thetahat))/8/omega1/a1hat;

Ixi  = 1e-4;
Ieta = 1e-4;

M  = 1000;
dt = 10^(-4);
N  = 2000000;
t  = (0:dt:N*dt)';

dist_opt_l = zeros(M,1);
phi        = zeros(M,N+1);
for i=1:M
    a1    = a1hat;
    T     = That;
    theta = thetahat;
    phi1  = phi1hatp;
    Xem = [a1 T theta phi1; zeros(N,4)];

    for j=1:N
        dW = sqrt(dt)*randn(1,2);
        Ia11   = (omega1^2)*(a1^2)*Ieta/4 + Ixi/(omega1^2);
        Iphi11 = 3*(omega1^2)*Ieta/4 + Ixi/((a1^2)*(omega1^2));
        Fnoise = [sqrt(Ia11)*dW(1); sqrt(Iphi11)*dW(2)];
        f = [-Gamma1*a1+(S*sin(Delta)-T*sin(theta))/2/omega1; ...
            a1-a1hat+That-T; ...
            a1-a1hat+2*T-2*That; ...
            (3*gamma*a1^3-4*S*cos(Delta)+4*T*cos(theta))/8/omega1/a1];
        a1    = a1    + dt*f(1) + Fnoise(1);
        T     = T     + dt*f(2);
        theta = theta + dt*f(3) + k*Fnoise(2);
        phi1  = phi1  + dt*f(4) + Fnoise(2);
        Xem(j+1,:) = [a1 T theta phi1];
    end
    phi(i,:) = Xem(:,4)';
    dist_opt_l(i) = Xem(end,4);
end

y5 = -0.01:0.0001:0.01;
mu5 = 0;
sigma5 = sqrt(var(dist_opt_l/sqrt(200)));
f5 = exp(-(y5-mu5).^2./(2*sigma5^2))./(sigma5*sqrt(2*pi));

figure(6)
histogram(dist_opt_l/sqrt(200)-mean(dist_opt_l/sqrt(200)), ...
    'LineWidth', 2, 'EdgeColor', [1 1 1], 'FaceColor',[1 0 0], ...
    'Normalization','pdf', 'NumBins',10)
hold on
plot(y5, f5, 'r--', 'LineWidth', 2)
xlim([-0.01,0.01])
set(gca, 'FontSize', 22)
ylabel('$w(\delta\phi(s)/\sqrt{s})$', 'Interpreter', 'Latex', ...
    'Fontsize', 24)
xlabel('$\delta\phi(s)/\sqrt{s}$', 'Interpreter', 'Latex', ...
    'Fontsize', 24)
