%% Read the s2p files
fns = 'forward.s2p';
[f(:),S11_real(:),S11_imag(:),S21_real(:),S21_imag(:),six(:),seven(:),eight(:),nine(:)] = textread(fns, '%f%f%f%f%f%f%f%f%f', 'headerlines', 1);
S11 = complex(S11_real, S11_imag);
S21 = complex(S21_real, S21_imag);

fns2 = 'reverse.s2p';
[f(:),S22_real(:),S22_imag(:),S12_real(:),S12_imag(:),six(:),seven(:),eight(:),nine(:)] = textread(fns2, '%f%f%f%f%f%f%f%f%f', 'headerlines', 1);
S22 = complex(S22_real, S22_imag);
S12 = complex(S12_real, S12_imag);
%% Plotting S11 Magnitude and Phase
S11_mag = mag2db(abs(S11));
S11_phase = rad2deg(angle(S11));

figure(1);
title('S11 Scattering Matrix')

yyaxis left
plot(f/1E6,S11_mag)
ylim([-20,5])
ylabel('S11 Magnitude (in dB)')

hold on
yyaxis right
plot(f/1E6,S11_phase)
ylim([-150,150])
ylabel('S11 Phase (in Degrees)')
xlabel('Frequency (MHz)');
legend('S11 Magnitude', 'S11 Phase');

xlim([0, 400]);
xticks(0:100:400)
%% Plotting S21 Magnitude and Phase
S21_mag = mag2db(abs(S21));
S21_phase = rad2deg(angle(S21));

figure(2);
title('S21 Scattering Matrix')

xlim([0, 400]);
yyaxis left
plot(f/1E6,S21_mag)
ylim([-20,5])
ylabel('S21 Magnitude (in dB)')
xlabel('Frequency (MHz)');

hold on
yyaxis right
ylim([-150,150])
plot(f/1E6,S21_phase)
ylabel('S21 Phase (in Degrees)')
legend('S21 Magnitude', 'S21 Phase');
%% Plotting S22 Magnitude and Phase
S22_mag = mag2db(abs(S22));
S22_phase = rad2deg(angle(S22));

figure(3);
title('S22 Scattering Matrix')

xlim([0, 400]);
yyaxis left
plot(f/1E6,S22_mag)
ylim([-20,5])
ylabel('S22 Magnitude (in dB)')
xlabel('Frequency (MHz)');

hold on
yyaxis right
plot(f/1E6,S22_phase)
ylim([-150,150])
ylabel('S22 Phase (in Degrees)')
legend('S22 Magnitude', 'S22 Phase');
%% Plotting S12 Magnitude and Phase
S12_mag = mag2db(abs(S12));
S12_phase = rad2deg(angle(S12));

figure(4);
title('S12 Scattering Matrix')

xlim([0, 400]);
yyaxis left
plot(f/1E6,S12_mag)
ylim([-20,5])
ylabel('S12 Magnitude (in dB)')
xlabel('Frequency (MHz)');

hold on
yyaxis right
ylim([-150,150])
plot(f/1E6,S12_phase)
ylabel('S12 Phase (in Degrees)')
legend('S12 Magnitude', 'S12 Phase');
%% Plotting Z11
Z11 = 50.*((((1+S11).*(1-S22))+(S12.*S21))./(((1-S11).*(1-S22))-(S12.*S21)));

figure(5);
hold on;

plot(f/1E6, real(Z11), '.r', 'Markersize', 5)
plot(f/1E6, imag(Z11), '.b', 'Markersize', 5)
xlabel('Frequency (MHz)');
ylabel('Z11 (Ohms)');
title('Z11 Impedance Matrix');
legend('Real(Z11)', 'Imag(Z11)');

xlim([0, 400]);
ylim([0, 100])
yticks(0:20:100)
xticks(0:100:400)
%% Plotting Z12
Z12 = 50.*((2.*S12)./((1-S11).*(1-S22)-(S12.*S21)));

figure(6);
hold on;

plot(f/1E6, real(Z12), '.r', 'Markersize', 5)
plot(f/1E6, imag(Z12), '.b', 'Markersize', 5)
xlabel('Frequency (MHz)');
ylabel('Z12 (Ohms)');
title('Z12 Impedance Matrix');
legend('Real(Z12)', 'Imag(Z12)');

xlim([0, 400]);
ylim([0, 100])
xticks(0:100:400)
yticks(0:20:100)
%% Plotting Z12
Z21 = 50.*((2.*S21)./((1-S11).*(1-S22)-(S12.*S21)));

figure(7);
hold on;

plot(f/1E6, real(Z21), '.r', 'Markersize', 5)
plot(f/1E6, imag(Z21), '.b', 'Markersize', 5)
xlabel('Frequency (MHz)');
ylabel('Z21 (Ohms)');
title('Z21 Impedance Matrix');
legend('Real(Z21)', 'Imag(Z21)');

xlim([0, 400]);
ylim([0, 100])
xticks(0:100:400)
yticks(0:20:100)
%% Plotting Z22
Z22 = 50.*((((1-S11).*(1+S22))+(S12.*S21))./(((1-S11).*(1-S22))-(S12.*S21)));

figure(8);

hold on;
plot(f/1E6, real(Z22), '.r', 'Markersize', 5)
plot(f/1E6, imag(Z22), '.b', 'Markersize', 5)
xlabel('Frequency (MHz)');
ylabel('Z22 (Ohms)');
title('Z22 Impedance Matrix');
legend('Real(Z22)', 'Imag(Z22)');

xlim([0, 400]);
ylim([0, 100])
xticks(0:100:400)
yticks(0:20:100)