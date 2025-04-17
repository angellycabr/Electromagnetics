% Read the s1p files
fns = 'shorted_cable.s1p';
[f(:),S11_real(:),S11_imag(:)] = textread(fns, '%f %f %f', 'headerlines', 1);
S11_short = complex(S11_real, S11_imag);

fns = 'opened_cable.s1p';
[f(:),S12_real(:),S12_imag(:)] = textread(fns, '%f %f %f', 'headerlines', 1);
S11_open = complex(S12_real, S12_imag);
%% For Z0
Zis = 50*(1+S11_short)./(1-S11_short);
Zio = 50*(1+S11_open)./(1-S11_open);
Z0 = sqrt(Zis .* Zio);

figure;
hold on;
plot(f/1E6, real(Z0), '.r', 'Markersize', 5)
plot(f/1E6, imag(Z0), '.b', 'Markersize', 5)
xlabel('Frequency (MHz)');
ylabel('Z0 (Ohms)');
title('Characteristic Impedance');
legend('Real(Zc)', 'Imag(Zc)');

ylim([-10, 70])
xlim([0,1000])
xticks(0:200:1000)
%% For VF
gamma = log((sqrt(Zio ./ Zis) + 1) ./ (sqrt(Zio ./ Zis) - 1));
beta = unwrap(imag(gamma)) ./ 4;
c0 = 3*1E8;

VF = (2*pi.*f) ./ (c0 .* beta);

figure(2)
plot(f/1E6, VF, '-')
xlabel('Frequency (MHz)');
ylabel('Velocity Factor up/c0');
title('Velocity Factor');

xlim([0,1500])
ylim([0.6,0.8])
yticks(0.6:0.05:0.8)
%xticks(0:500:1500)
%% For Alpha
alpha = unwrap(real(gamma)) ./ 4; 

figure
plot(f/10^6, alpha, '-')
xlabel('Frequency (MHz)');
ylabel('Alpha (Np/m)');
title('Attenuation Constant');

ylim([-0.1,0.5])
xlim([0,1000])
xticks(0:200:1000)
%% Applying Polyfit for Alpha Plots
p = polyfit(f/10^6, alpha, 2);

% Create a range of frequencies for plotting the fitted curve
f_range = linspace(0, 6000, 10000); % Adjust the number of points as needed

% Evaluate the polynomial at the specified frequency range
alpha_fit = polyval(p, f_range);

% Plot the original data and the fitted curve
figure
plot(f/10^6, alpha, '-', f_range, alpha_fit, '--')
xlabel('Frequency (MHz)');
ylabel('Alpha (Np/m)');
title('Attenuation Constant');
legend('Original Data', 'Fitted Curve');
ylim([-0.1, 0.5])
xlim([0, 1000])
xticks(0:200:1000)
%% Applying Polyfit for Z0 Plots
p_real = polyfit(f/1E6, real(Z0), 2);
p_imag = polyfit(f/1E6, imag(Z0), 2);

f_range = linspace(0, 6000, 10000); 
Z0_real_fit = polyval(p_real, f_range);
Z0_imag_fit = polyval(p_imag, f_range);

% Plot the original data and the fitted curves
figure;
hold on;
plot(f/1E6, real(Z0), '.r', 'Markersize', 3)
plot(f/1E6, imag(Z0), '.b', 'Markersize', 3)
plot(f_range, Z0_real_fit, '-r')
plot(f_range, Z0_imag_fit, '-b')
xlabel('Frequency (MHz)');
ylabel('Z0 (Ohms)');
title('Characteristic Impedance');
legend('Real(Z0)', 'Imag(Z0)', 'Fitted Real(Z0)', 'Fitted Imag(Z0)');
ylim([-10, 70])
xlim([0, 1000])
xticks(0:200:1000)
