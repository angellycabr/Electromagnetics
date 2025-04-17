clear;
clc;
close all;
%% 6.2
fns = 'Lab 6/unbalanced2.s1p';
[f(:),S11_real(:),S11_imag(:)] = textread(fns, '%f %f %f', 'headerlines', 1);
S11 = complex(S11_real, S11_imag);

w = 2*pi.*f;
Cs = 1.09*1E-12; %Inductance of capacitor
Ls = 3.41*1E-9; %Inductance of inductance

Zi = 50*(1+S11)./(1-S11);
S11 = (Zi - 50) ./ (Zi + 50); 

S11_mag = mag2db(abs(S11));
S11_phase = (180/pi) .* angle(S11);

windowSize = 30;
S11_mag_smoothed = smoothdata(S11_mag, 'movmean', windowSize);
S11_phase_smoothed = smoothdata(S11_phase, 'movmean', windowSize);


%S11_phase = rad2deg(angle(S11));
%rad2deg(angle(S11));

figure;
title('S11 Unbalanced Dipole Magntiude and Phase');
yyaxis left
plot(f/1E6, S11_mag_smoothed);
xlim([50, 1000])
%xticks(0.05:50:650)
ylim([-20, 5])
%yticks(-20:10:5)
xlabel('Frequency (MHz)')
ylabel('|S11| (in dB)');
    
hold on
yyaxis right
plot(f/1E6, S11_phase_smoothed);
ylim([-150, 150]);
%yticks(-150:50:150)
ylabel('<S11 (in deg)');

windowSize = 25;
Zi_smoothed = smoothdata(Zi, 'movmean', windowSize);

figure;
title('Unbalanced Dipole Input Impedance Zi=R+jX');
yyaxis left
plot(f/1E6, real(Zi_smoothed));
xlim([500, 700])
xticks(500:50:700)
ylim([-100, 300])
yticks(-100:100:300)
xlabel('Frequency (MHz)')
ylabel('R (in ohms)');
    
hold on
yyaxis right
plot(f/1E6, imag(Zi_smoothed));
ylim([-1000, 1000]);
yticks(-1000:500:1000)
ylabel('X (in ohms)');
%% 6.4
clear;
clc;
close all;

fns = 'Lab 6/balanced.s1p';
[f(:),S11_real(:),S11_imag(:)] = textread(fns, '%f %f %f', 'headerlines', 1);
S11 = complex(S11_real, S11_imag);

S11_mag = mag2db(abs(S11));
S11_phase = (180/pi) .* angle(S11);

windowSize = 100;
S11_mag_smoothed = smoothdata(S11_mag, 'movmean', windowSize);
S11_phase_smoothed = smoothdata(S11_phase, 'movmean', windowSize);

%S11_phase = rad2deg(angle(S11));
%rad2deg(angle(S11));

figure;
title('S11 Balanced Dipole Magntiude and Phase');
yyaxis left
plot(f/1E6, S11_mag_smoothed);
xlim([50, 1000])
%xticks(0.05:50:650)
ylim([-20, 5])
%yticks(-20:10:5)
xlabel('Frequency (MHz)')
ylabel('|S11| (in dB)');
    
hold on
yyaxis right
plot(f/1E6, S11_phase_smoothed);
ylim([-150, 150]);
%yticks(-150:50:150)
ylabel('<S11 (in deg)');
%%
clear;
clc;
close all;

fns = 'Lab 6/balanced2.s2p';
[f(:),S11_real(:),S11_imag(:),S21_real(:),S21_imag(:),six(:),seven(:),eight(:),nine(:)] = textread(fns, '%f%f%f%f%f%f%f%f%f', 'headerlines', 1);
S11 = complex(S11_real, S11_imag);
S21 = complex(S21_real, S21_imag);

S11_mag = mag2db(abs(S11));
S11_phase = (180/pi) .* angle(S11);

S21_mag = mag2db(abs(S21));
S21_phase = (180/pi) .* angle(S21);

windowSize = 50;
S11_mag_smoothed = smoothdata(S11_mag, 'movmean', windowSize);
S11_phase_smoothed = smoothdata(S11_phase, 'movmean', windowSize);

windowSize = 100;
S21_mag_smoothed = smoothdata(S21_mag, 'movmean', windowSize);
%S21_phase_smoothed = smoothdata(S21_phase, 'movmean', windowSize);

S21_phase = rad2deg(angle(S21));
%rad2deg(angle(S11));

figure;
title('S11 Balanced Dipole Magntiude and Phase');
yyaxis left
plot(f/1E6, S11_mag_smoothed);
xlim([50, 1000])
%xticks(0.05:50:650)
ylim([-20, 5])
%yticks(-20:10:5)
xlabel('Frequency (MHz)')
ylabel('|S11| (in dB)');
    
hold on
yyaxis right
plot(f/1E6, S11_phase_smoothed);
ylim([-150, 150]);
%yticks(-150:50:150)
ylabel('<S11 (in deg)');

figure;
title('S21 Balanced Dipole Magntiude and Phase');
yyaxis left
plot(f/1E6, S21_mag_smoothed);
xlim([50, 1000])
%xticks(0.05:50:650)
%ylim([-20, 5])
%yticks(-20:10:5)
xlabel('Frequency (MHz)')
ylabel('|S21| (in dB)');
    
hold on
yyaxis right
plot(f/1E6, S21_phase);
%ylim([-150, 150]);
%yticks(-150:50:150)
ylabel('<S21 (in deg)');

w = 2*pi.*f;
Cs = 1.09*1E-12; %Inductance of capacitor
Ls = 3.41*1E-9; %Inductance of inductance

Zi = 50*(1+S11)./(1-S11);

windowSize = 20;
Zi_smoothed = smoothdata(Zi, 'movmean', windowSize);

figure;
title('Balanced Dipole Input Impedance Zi=R+jX');
yyaxis left
plot(f/1E6, real(Zi_smoothed));
xlim([200, 1000])
xticks(200:200:1000)
ylim([0, 400])
yticks(0:100:400)
xlabel('Frequency (MHz)')
ylabel('R (in ohms)');
    
hold on
yyaxis right
plot(f/1E6, imag(Zi_smoothed));
ylim([-1000, 1000]);
yticks(-1000:500:1000)
ylabel('X (in ohms)');
%%
S22 = S11;
S12 = S21;

Z11 = 50.*((((1+S11).*(1-S22))+(S12.*S21))./(((1-S11).*(1-S22))-(S12.*S21)));
Z21 = 50.*((2.*S21)./((1-S11).*(1-S22)-(S12.*S21)));

Zi = 2.*(Z11 - Z21);

figure;
title('Balanced Dipole Input Impedance Zi=R+jX');
yyaxis left
plot(f/1E6, real(Zi));
xlim([500, 700])
xticks(500:50:700)
ylim([50, 250])
yticks(50:25:250)
xlabel('Frequency (MHz)')
ylabel('R (in ohms)');
    
hold on
yyaxis right
plot(f/1E6, imag(Zi));
ylim([-250, 50]);
yticks(-250:50:50)
ylabel('X (in ohms)');
%%
