fns = 'reflect.s2p';
[f(:),S11_real(:),S11_imag(:),S21_real(:),S21_imag(:),six(:),seven(:),eight(:),nine(:)] = textread(fns, '%f%f%f%f%f%f%f%f%f', 'headerlines', 1);
S11 = complex(S11_real, S11_imag);
S21 = complex(S21_real, S21_imag);
%% Plotting Gamma
S11_mag = mag2db(abs(S11));
S11_phase = rad2deg(angle(S11));

figure(1);
plot(f/1E6,S11_mag)
title('\Gamma2 Measurement')
ylabel('\Gamma2 Magnitude (in dB)')
xlabel('Frequency (MHz)');
%% Reading for Port 2
fns = 'reflectport2.s2p';
[f(:),S11_real(:),S11_imag(:),S21_real(:),S21_imag(:),six(:),seven(:),eight(:),nine(:)] = textread(fns, '%f%f%f%f%f%f%f%f%f', 'headerlines', 1);
S11 = complex(S11_real, S11_imag);
S21 = complex(S21_real, S21_imag);
%% Plotting Gamma
S11_mag = mag2db(abs(S11));
S11_phase = rad2deg(angle(S11));

figure(2);
plot(f/1E6,S11_mag)
title('\Gamma2 Measurement of Port 2')
ylabel('\Gamma2 Magnitude (in dB)')
xlabel('Frequency (MHz)');