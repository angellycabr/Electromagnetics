%% 4.5 and 4.6 Recalculate for Monopole
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

    Z11(i) = ((-1 / (sin(beta .* h)).^2) .* trapz(z,In1));
    Z21(i) = ((-1 / (sin(beta .* h)).^2) .* trapz(z,In2));

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
%% 4.5
Grt = 3.3;  % Monopole gain in linear scale
Gm = zeros(1, length(freq));

for i = 1:length(freq)
    lambda = (3 * 1E8) ./ freq(i);
    % Adjust S11 and S21 for monopole
    Gm(i) = (1 - abs(S11c(i)).^2) * Grt * sqrt((lambda / (4*pi*r2)).^2);
end

S21_magf = mag2db(abs(Gm));

figure;
title('S21 of Monopoles (Calculated with Modified Friis)');
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
%% 4.6
lambda = (3 * 1E8) ./ (600 * 1E6);
Gmes = sqrt(0.015 / ((lambda / (4*pi*r2)).^2));

Gm = zeros(1, length(freq));

for i = 1:length(freq)
    lambda = (3 * 1E8) ./ freq(i);
    % Adjust S11 and S21 for monopole
    Gm(i) = (1 - abs(S11c(i)).^2) * Gmes * sqrt((lambda / (4*pi*r2)).^2);
end

S21_magf = mag2db(abs(Gm));

figure;
title('S21 of Monopoles (Calculated Gm)');
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