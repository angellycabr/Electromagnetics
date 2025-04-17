clear;
clc;
close all;
%% 4.10
r = 0.0002;
h = 0.120; %250 / 2000;

% Specify z range
npoints = 500;
zstrt = 0;
zstop = h;
zstep = (zstop - zstrt) / (npoints - 1);
z = -zstop:zstep:zstop;

% Specify frequencies
nfreq = 1010;
fstrt = 50*1E3;
fstop = 1*1E9;
fstep = (fstop - fstrt) / (nfreq - 1);
freq = fstrt:fstep:fstop;
lambda = (3 * 1E8) ./ freq;

% Define parameters
eta = 120 * pi;
Im = 1;

R0 = sqrt(r^2 + z.^2);
R1 = sqrt(r^2 + (z - h).^2);
R2 = sqrt(r^2 + (z + h).^2);

Zi = zeros(1, length(freq));
S11 = zeros(1, length(freq));
    
for i = 1:length(freq)
    beta = (2 * pi) ./ lambda(i);
    E_z = ((1j * eta) / (4 * pi)) .* (((2 * exp(-1j .* beta * R0) * cos(beta * h)) ./ R0) - (exp(-1j .* beta * R1) ./ R1) - (exp(-1j .* beta * R2) ./ R2)); 
    H = (-1j * Im) / (4 * pi) * (1 / r) .* (2 * exp(-1j .* beta * R0) .* cos(beta * h) - exp(-1j .* beta * R1) - exp(-1j .* beta * R2));
        
    %In = E_z .* conj(H);
    %Zi(i) = ((-2 * pi * r) / (sin(beta * h)).^2) * trapz(z, In); %Eq. 25

    In = E_z .* sin(beta .* (h - abs(z)));
    Zi(i) = (-1 / (sin(beta .* h)).^2) .* trapz(z,In); %Eq. 30

    S11(i) = (Zi(i) - 50) / (Zi(i) + 50);
end

figure;
title(['S11 Antenna Magntiude and Phase (2h=' num2str(h*1000) ' mm)']);
yyaxis left
plot(freq/1E6, mag2db(abs(S11)));
%xlim([50, 100])
%xticks(0.05:50:650)
ylim([-20, 5])
%yticks(-20:10:5)
xlabel('Frequency (MHz)')
ylabel('|S11| (in dB)');
    
hold on
yyaxis right
plot(freq/1E6, rad2deg(angle(S11)));
ylim([-150, 150]);
%yticks(-150:50:150)
ylabel('<S11 (in deg)');
%% 4.11
fns = 'Lab 6/antenna2.s1p';
[f(:),S11_real(:),S11_imag(:)] = textread(fns, '%f %f %f', 'headerlines', 1);
S11 = complex(S11_real, S11_imag);

S11_mag = mag2db(abs(S11));
S11_phase = (180/pi) .* angle(S11);

windowSize = 12;
S11_mag_smoothed = smoothdata(S11_mag, 'movmean', windowSize);
S11_phase_smoothed = smoothdata(S11_phase, 'movmean', windowSize);

%S11_phase = rad2deg(angle(S11));
%rad2deg(angle(S11));

figure;
title('S11 Antenna Measured Magntiude and Phase');
yyaxis left
plot(f/1E6, S11_mag_smoothed);
xlim([50, 1000])
%xticks(0.05:50:650)
ylim([-20, 5])
%yticks(-20:10:5)
xlabel('Frequency (MHz)')
ylabel('|S11| (in deg)');
    
hold on
yyaxis right
plot(f/1E6, S11_phase_smoothed);
ylim([-150, 150]);
%yticks(-150:50:150)
ylabel('<S11 (in deg)');
%% 4.11 Calculated S11 of Antenna
S11_calc = (0.5*Zi-50)./(0.5*Zi+50);
Ze = (50 * (1 + S11)) ./ (1 - S11);
Zec = (50 * (1 + S11_calc)) ./ (1 - S11_calc);

S11_magc = mag2db(abs(S11_calc));
S11_phasec = rad2deg(angle(S11_calc));

figure;
title('S11 Monopole Magntiude and Phase');
yyaxis left
plot(freq/1E6, S11_magc);
xlim([0, 1000])
xticks(0:200:1000)
ylim([-20, 5])
%yticks(-20:10:5)
xlabel('Frequency (MHz)')
ylabel('|S11| (in dB)');
    
hold on
yyaxis right
plot(freq/1E6, S11_phasec);
ylim([-150, 150]);
%yticks(-150:50:150)
ylabel('<S11 (in dB)');
%% 5.2
w = 2*pi.*freq;
Cs = 1.09*1E-12; %Inductance of capacitor
Ls = 3.41*1E-9; %Inductance of inductance

Zi = (1 ./ ((1./Ze)-1j.*w.*Cs))-1j.*w.*Ls;
S11_new = (Zi - 50) ./ (Zi + 50);
S11_mag_new = mag2db(abs(S11_new));
S11_phase_new = rad2deg(angle(S11_new));

figure;
title('S11 Monopole Magntiude and Phase');
yyaxis left
plot(f/1E6, S11_mag_new);
hold on
plot(freq/1E6, S11_magc);
xlim([0, 1000])
xticks(0:200:1000)
ylim([-20, 5])
%yticks(-20:10:5)
xlabel('Frequency (MHz)')
ylabel('|S11| (in dB)');
    
hold on
yyaxis right
plot(f/1E6, S11_phase_new);
plot(freq/1E6, S11_phasec);
ylim([-150, 150]);
%yticks(-150:50:150)
ylabel('<S11 (in deg)');
legend("Measured", "Calculated", "Measured", "Calculated")

figure;
title('Dipoles Input Impedance, Zi=R+jX (Predicted)');
yyaxis left
plot(f/ 1E6, real(Zi));
xlim([500, 700])
xticks(500:50:700)
ylim([-100, 300])
yticks(-100:100:300)
xlabel('Frequency (MHz)')
ylabel('R (in ohms)');
    
hold on
yyaxis right
plot(f/ 1E6, imag(Zi));
ylim([-1000, 1000]);
yticks(-1000:500:1000)
ylabel('X (in ohms)');