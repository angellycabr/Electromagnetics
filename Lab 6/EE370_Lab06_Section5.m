clear;
clc;
close all;
%% 5.2 
fns = 'Lab 6/antenna2.s1p';
[f(:),S11_real(:),S11_imag(:)] = textread(fns, '%f %f %f', 'headerlines', 1);
S11 = complex(S11_real, S11_imag);

Ze = 50.*((1 + S11) ./ (1- S11));

w = 2*pi.*f;
Cs = 1.09*1E-12; %Inductance of capacitor
Ls = 3.41*1E-9; %Inductance of inductance
%Zi = 1./((1./Ze) - (1j.*w.*Cs)) - (1j.*w.*Ls);
%Zi = Ze;
Zi = (1 ./ ((1./Ze)-1j.*w.*Cs))-1j.*w.*Ls;

S11_new = (Zi - 50) ./ (Zi + 50);
S11_mag = mag2db(abs(S11_new));
S11_phase = rad2deg(angle(S11_new));

figure;
title('S11 Monopole Magntiude and Phase');
yyaxis left
plot(f/1E6, S11_mag);
xlim([0, 1000])
xticks(0:200:1000)
ylim([-20, 5])
%yticks(-20:10:5)
xlabel('Frequency (MHz)')
ylabel('|S11| (in dB)');
    
hold on
yyaxis right
plot(f/1E6, S11_phase);
ylim([-150, 150]);
%yticks(-150:50:150)
ylabel('<S11 (in dB)');

figure;
title('S11 Monopole Magntiude and Phase');
yyaxis left
plot(f/1E6, real(Ze));
%xlim([0, 1000])
%xticks(0:200:1000)
ylim([-20, 5])
%yticks(-20:10:5)
xlabel('Frequency (MHz)')
ylabel('|S11| (in dB)');
    
hold on
yyaxis right
plot(f/1E6, imag(Ze));
ylim([-100, 100]);
yticks(-100:50:100)
ylabel('<S11 (in dB)');
%% 5.6
clear;
clc;
close all;

fns = 'Lab 6/90degantenna.s1p';
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
title('S11 90 Degree Antenna Measured Magntiude and Phase');
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
%% 5.6
clear;
clc;
close all;

fns = 'Lab 6/135degantenna.s1p';
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
title('S11 135 Degree Antenna Measured Magntiude and Phase');
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