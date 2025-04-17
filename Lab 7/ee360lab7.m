clc;
clear;
close all;
%% 2.1
r = 0.0002;
h = 0.120;
r2 = 0.5;

% Specify z range
npoints = 400;
zstrt = 0;
zstop = h;
zstep = (zstop - zstrt) / (npoints - 1);
z = -zstop:zstep:zstop;

% Specify frequencies
nfreq = 400;
fstrt = 50*1E3;
fstop = 1*1E9;
fstep = (fstop - fstrt) / (nfreq - 1);
freq = fstrt:fstep:fstop;
lambda = (3 * 1E8) ./ freq;
beta = (2 * pi) ./ lambda;

% Define parameters
eta = 120 * pi;
Im = 1;

R0 = sqrt(r^2 + z.^2);
R1 = sqrt(r^2 + (z - h).^2);
R2 = sqrt(r^2 + (z + h).^2);

R02 = sqrt(r2^2 + z.^2);
R12 = sqrt(r2^2 + (z - h).^2);
R22 = sqrt(r2^2 + (z + h).^2);

Z11 = zeros(1, length(freq));
Z21 = zeros(1, length(freq));
S11c = zeros(1, length(freq));
S21c = zeros(1, length(freq));
    
for i = 1:length(freq)
    beta = (2 * pi) ./ lambda(i);
    E_z1 = ((1j * eta) / (4 * pi)) .* (((2 * exp(-1j .* beta * R0) * cos(beta * h)) ./ R0) - (exp(-1j .* beta * R1) ./ R1) - (exp(-1j .* beta * R2) ./ R2)); 
    E_z2 = ((1j * eta) / (4 * pi)) .* (((2 * exp(-1j .* beta * R02) * cos(beta * h)) ./ R02) - (exp(-1j .* beta * R12) ./ R12) - (exp(-1j .* beta * R22) ./ R22)); 

    In1 = E_z1 .* sin(beta .* (h - abs(z)));
    In2 = E_z2 .* sin(beta .* (h - abs(z)));

    Z11(i) = ((-1 / (sin(beta .* h)).^2) .* trapz(z,In1)) ./2;
    Z21(i) = ((-1 / (sin(beta .* h)).^2) .* trapz(z,In2)) ./2;

    S11c(i) = (Z11(i) - 50) / (Z11(i) + 50);
    S21c(i) = (2 * Z21(i) * 50) ./ ((Z11(i) + 50).^2 - Z21(i).^2);
end

figure;
title('S11 Dipole Magntiude and Phase (Calculated)');
yyaxis left
plot(freq/1E6, mag2db(abs(S11c)));
xlim([0, 1000])
xticks(0:200:1000)
ylim([-30, 5])
yticks(-30:5:5)
xlabel('Frequency (MHz)')
ylabel('|S11| (in dB)');
    
hold on
yyaxis right
plot(freq/1E6, rad2deg(angle(S11c)));
ylim([-150, 150]);
yticks(-150:50:150)
ylabel('<S11 (in deg)');

figure;
title('S21 Dipole Magntiude and Phase (Calculated)');
yyaxis left
plot(freq/1E6, mag2db(abs(S21c)));
xlim([0, 1000])
xticks(0:200:1000)
ylim([-120, 0])
yticks(-120:20:0)
xlabel('Frequency (MHz)')
ylabel('|S21| (in dB)');
    
hold on
yyaxis right
plot(freq/1E6, rad2deg(angle(S21c)));
ylim([-200, 200]);
%yticks(-150:50:150)
ylabel('<S21 (in deg)');
%% 2.3
%clear;
close all;

fns = 'Lab 7/communicate.s2p';
[f(:),S11_real(:),S11_imag(:),S21_real(:),S21_imag(:),six(:),seven(:),eight(:),nine(:)] = textread(fns, '%f%f%f%f%f%f%f%f%f', 'headerlines', 1);
S11 = complex(S11_real, S11_imag);
S21 = complex(S21_real, S21_imag);

S11_mag = mag2db(abs(S11));
S11_phase = rad2deg(angle(S11));

S21_mag = mag2db(abs(S21));
S21_phase = rad2deg(angle(S21));

figure;
title('S11 of Dipoles');
yyaxis left
plot(f/1E6, S11_mag);
hold on
plot(freq/1E6, mag2db(abs(S11c)));
xlim([0, 1000])
xticks(0:200:1000)
ylim([-30, 5])
yticks(-30:5:5)
xlabel('Frequency (MHz)')
ylabel('|S11| (in dB)');
    
hold on
yyaxis right
plot(f/1E6, S11_phase);
plot(freq/1E6, rad2deg(angle(S11c)));
ylim([-150, 150]);
yticks(-150:50:150)
ylabel('<S11 (in deg)');
legend("Measured", "Calculated", "Measured", "Calculated")

figure;
title('S21 of Dipoles');
yyaxis left
plot(f/1E6, S21_mag);
hold on;
plot(freq/1E6, mag2db(abs(S21c)));
xlim([0, 1000])
xticks(0:200:1000)
ylim([-80, 0])
yticks(-80:20:0)
xlabel('Frequency (MHz)')
ylabel('|S21| (in dB)');
    
hold on
yyaxis right
plot(f/1E6, S21_phase);
plot(freq/1E6, rad2deg(angle(S21c)));
ylim([-150, 150]);
yticks(-150:50:150)
ylabel('<S21 (in deg)');
legend("Measured", "Calculated", "Measured", "Calculated")
%% 3.3
% Pre-defined Parameters
Cs = 1.09*1E-12;
Ls = 3.41*1E-9;

S11_em = zeros(1, length(f));
S21_em = zeros(1, length(f));

for i=1:length(f)
    w = 2*pi*f(i);
    XL = 1i*w*Ls;
    BC = 1i*w*Cs;

    T11 = ((1 + S11(i)) .* (1 - S11(i)) + S21(i).^2) ./ (2 * S21(i));
    T12 = 50*(((1 + S11(i)).^2 - S21(i).^2) ./ (2 * S21(i)));
    T21 = 1/50*(((1 - S11(i)) .* (1 - S11(i)) - S21(i).^2) ./ (2 * S21(i)));
    T22 = ((1 - S11(i)) .* (1 + S11(i)) + S21(i).^2) ./ (2 * S21(i));

    Te = [T11, T12; T21, T22];
    T1 = [1, XL;BC, 1+(XL.*BC)];
    T2 = [1+(XL .* BC) XL; BC 1];
    T = inv(T1) * Te * inv(T2);

    S11_em(i) = (T(1,1) + T(1,2)/50 - T(2,1)*50 - T(2,2)) / (T(1,1) + T(1,2)/50 + T(2,1)*50 + T(2,2));
    S21_em(i) = 2/(T(1,1) + T(1,2)/50 + T(2,1)*50 + T(2,2));
end

figure;
title('S11 Dipole Magntiude and Phase (Measured)');
yyaxis left
plot(f/1E6, mag2db(abs(S11_em)));
hold on;
plot(freq/1E6, mag2db(abs(S11c)));
xlim([0, 1000])
xticks(0:200:1000)
ylim([-30, 5])
yticks(-30:5:5)
xlabel('Frequency (MHz)')
ylabel('|S11| (in dB)');
    
hold on
yyaxis right
plot(f/1E6, rad2deg(angle(S11_em)));
plot(freq/1E6, rad2deg(angle(S11c)));
ylim([-150, 150]);
yticks(-150:50:150)
ylabel('<S11 (in deg)');
legend("Measured", "Calculated", "Measured", "Calculated")

figure;
title('S21 Dipole Magntiude and Phase (Measured)');
yyaxis left
plot(f/1E6, mag2db(abs(S21_em)));
hold on;
plot(freq/1E6, mag2db(abs(S21c)));
xlim([0, 1000])
xticks(0:200:1000)
ylim([-120, 0])
yticks(-120:20:0)
xlabel('Frequency (MHz)')
ylabel('|S21| (in dB)');
    
hold on
yyaxis right
plot(f/1E6, rad2deg(angle(S21_em)));
plot(freq/1E6, rad2deg(angle(S21c)));
ylim([-150, 150]);
%yticks(-150:50:150)
ylabel('<S21 (in deg)');
legend("Measured", "Calculated", "Measured", "Calculated")
%% 4.3
Gr = 5.2; 
%lambda = (3 * 1E8) ./ (600 * 1E6);
%Grt = sqrt(0.015 / ((lambda / (4*pi*r2)).^2));

G = zeros(1, length(freq));

for i = 1:length(freq)
    lambda = (3 * 1E8) ./ freq(i);
    G(i) = 1.64 * (lambda / (4*pi*r2));
end

S21_magf = mag2db(abs(G));

figure;
title('S21 of Dipoles (Calculated with Friis)');
yyaxis left
plot(freq/1E6, mag2db(abs(S21c)));
hold on;
plot(freq/1E6, S21_magf);
xlim([0, 1000])
xticks(0:200:1000)
ylim([-80, 0])
yticks(-80:20:0)
xlabel('Frequency (MHz)')
ylabel('|S21| (in dB)');
    
hold on
yyaxis right
plot(freq/1E6, rad2deg(angle(S21c)));
ylim([-150, 150]);
yticks(-150:50:150)
ylabel('<S21 (in deg)');
legend("Original Calculation", "Calculated w/ Friis")
%% 4.4
%lambda = (3 * 1E8) ./ (600 * 1E6);
%Grt = sqrt(0.015 / ((lambda / (4*pi*r2)).^2));

G = zeros(1, length(freq));

for i = 1:length(freq)
    lambda = (3 * 1E8) ./ freq(i);
    G(i) = (1 - abs(S11c(i)).^2) * 1.64 * sqrt((lambda / (4*pi*r2)).^2);
end

S21_magf = mag2db(abs(G));

figure;
title('S21 of Dipoles (Calculated with Modified Friis)');
yyaxis left
plot(freq/1E6, mag2db(abs(S21c)));
hold on;
plot(freq/1E6, S21_magf);
xlim([0, 1000])
xticks(0:200:1000)
ylim([-80, 0])
yticks(-80:20:0)
xlabel('Frequency (MHz)')
ylabel('|S21| (in dB)');
    
hold on
yyaxis right
plot(freq/1E6, rad2deg(angle(S21c)));
ylim([-150, 150]);
yticks(-150:50:150)
ylabel('<S21 (in deg)');
legend("Original Calculation", "Calculated w/ Friis")
%% 4.5
Grt = 3.3;  % Monopole gain in linear scale
Gm = zeros(1, length(freq));

for i = 1:length(freq)
    lambda = (3 * 1E8) ./ freq(i);
    % Adjust S11 and S21 for monopole
    S11_monopole = S11c(i) / 2;
    Gm(i) = (1 - abs(S11_monopole).^2) * Grt * sqrt((lambda / (4*pi*r2)).^2);
end

S21_magf = mag2db(abs(Gm));

figure;
title('S21 of Monopoles (Calculated with Modified Friis)');
yyaxis left
plot(freq/1E6, mag2db(abs(S21c ./ 2)));
hold on;
plot(freq/1E6, S21_magf);
xlim([0, 1000])
xticks(0:200:1000)
ylim([-80, 0])
yticks(-80:20:0)
xlabel('Frequency (MHz)')
ylabel('|S21| (in dB)');

hold on
yyaxis right
plot(freq/1E6, rad2deg(angle(S21c ./ 2)));
ylim([-150, 150]);
yticks(-150:50:150)
ylabel('<S21 (in deg)');
legend("Original Calculation", "Calculated w/ Friis")
%% 4.6
% Measure Gm of Monopole
lambda = (3 * 1E8) ./ (600 * 1E6);
%Gmes = sqrt((0.15) * ((4 * pi * r2 / lambda)^2));
Gmes = sqrt(0.015 / ((lambda / (4*pi*r2)).^2));

Gm = zeros(1, length(freq));

for i = 1:length(freq)
    lambda = (3 * 1E8) ./ freq(i);
    % Adjust S11 and S21 for monopole
    S11_monopole = S11c(i) / 2;
    Gm(i) = (1 - abs(S11c(i)).^2) * 1.53 * sqrt((lambda / (4*pi*r2)).^2);
end

S21_magf = mag2db(abs(Gm));

figure;
title('S21 of Monopoles (Calculated with Modified Friis)');
yyaxis left
plot(freq/1E6, mag2db(abs(S21 ./ 2)));
hold on;
plot(freq/1E6, S21_magf);
xlim([0, 1000])
xticks(0:200:1000)
ylim([-80, 0])
yticks(-80:20:0)
xlabel('Frequency (MHz)')
ylabel('|S21| (in dB)');

hold on
yyaxis right
plot(freq/1E6, rad2deg(angle(S21 ./ 2)));
ylim([-150, 150]);
yticks(-150:50:150)
ylabel('<S21 (in deg)');
legend("Original Calculation", "Calculated w/ Friis")
