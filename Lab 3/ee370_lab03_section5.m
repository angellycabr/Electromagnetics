%% Read the s2p files
fns = '2000mm.s2p';
[f(:),S11_real(:),S11_imag(:),S21_real(:),S21_imag(:),six(:),seven(:),eight(:),nine(:)] = textread(fns, '%f%f%f%f%f%f%f%f%f', 'headerlines', 1);
S11 = complex(S11_real, S11_imag);
S21 = complex(S21_real, S21_imag);
%% Plotting S21 for 2000 mm
S21_mag = mag2db(abs(S21));
S21_phase = rad2deg(unwrap(angle(S21)));

figure(2);
title('S21 Scattering Matrix')

xlim([0, 1500]);
yyaxis left
plot(f/1E6,S21_mag)
ylim([-20,20])
ylabel('S21 Magnitude (in dB)')
xlabel('Frequency (MHz)');

hold on
yyaxis right
plot(f/1E6,S21_phase)
ylabel('S21 Phase (in Degrees)')
legend('S21 Magnitude', 'S21 Phase');
%% For Attenuation Alpha for 2000 mm
alpha = -20*log10(abs(S21)) ./ 4;

figure (1);
plot(f/10^6, alpha, '-')
xlabel('Frequency (MHz)');
ylabel('Alpha (in dB)');
title('Attenuation Constant');
%% For Velocity Factor for 2000 mm
c0 = 3*1E8;
VF = -1 .* (4*pi.*f) ./ (c0.*unwrap(angle(S21)));

figure(3)
plot(f/1E6, VF, '-')
xlabel('Frequency (MHz)');
ylabel('Velocity Factor up/c0');
title('Velocity Factor');

xlim([0,1500])
ylim([0.7,0.9])
yticks(0.7:0.05:0.9)
%% Repeat for 200 mm
fns2 = '200mm_section_5.s2p';
[f(:),S11_real(:),S11_imag(:),S21_real(:),S21_imag(:),six(:),seven(:),eight(:),nine(:)] = textread(fns2, '%f%f%f%f%f%f%f%f%f', 'headerlines', 1);
S11_200mm = complex(S11_real, S11_imag);
S21_200mm = complex(S21_real, S21_imag);
%% Plotting S21 for 200 mm
S21_200mm_mag = mag2db(abs(S21_200mm));
S21_200mm_phase = unwrap(angle(S21_200mm));

figure(3);
title('S21 Scattering Matrix for 200mm Cable')

xlim([0,1500])
yyaxis left
plot(f/1E6,S21_200mm_mag)
ylim([-20,20])
ylabel('S21 Magnitude (in dB)')
xlabel('Frequency (MHz)');

hold on
yyaxis right
plot(f/1E6,S21_200mm_phase)
ylabel('S21 Phase (in Degrees)')
legend('S21 Magnitude', 'S21 Phase');
%% For Attenuation Alpha of 200 mm
alpha = -20*log10(abs(S21_200mm)) ./ (2*0.21);

figure (4);
plot(f/10^6, alpha, '-')
xlabel('Frequency (MHz)');
ylabel('Alpha (in dB)');
title('Attenuation Constant for 200mm Cable');
%% For VF of 200 mm
c0 = 3*1E8;
VF = -1 .* (2*0.21*pi.*f) ./ (c0.*unwrap(angle(S21_200mm)));

figure(5)
plot(f/1E6, VF, '-')
xlabel('Frequency (MHz)');
ylabel('Velocity Factor up/c0');
title('Velocity Factor');

xlim([0,1500])
ylim([0.5,0.8])
yticks(0.5:0.05:0.8)