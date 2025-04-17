%% Read the s1p files for 210 mm cable
fns = '210mm_shorted_cable.s1p';
[f(:),S11_real(:),S11_imag(:)] = textread(fns, '%f %f %f', 'headerlines', 1);
S11_short = complex(S11_real, S11_imag);

fns2 = '210mm_opened_cable.s1p';
[f(:),S12_real(:),S12_imag(:)] = textread(fns2, '%f %f %f', 'headerlines', 1);
S11_open = complex(S12_real, S12_imag);
%% For Z0
Zis = 50*(1+S11_short)./(1-S11_short);
Zio = 50*(1+S11_open)./(1-S11_open);
Z0 = sqrt(Zis .* Zio);

figure (1);
hold on;
plot(f/1E6, real(Z0), '.r', 'Markersize', 5)
plot(f/1E6, imag(Z0), '.b', 'Markersize', 5)
xlabel('Frequency (MHz)');
ylabel('Z0 (Ohms)');
title('Characteristic Impedance');
legend('Real(Zc)', 'Imag(Zc)');

%ylim([-10, 70])
%xlim([0,1000])
%xticks(0:200:1000)
%% For VF
gamma = log((sqrt(Zio ./ Zis) + 1) ./ (sqrt(Zio ./ Zis) - 1));
beta = unwrap(imag(gamma))./(2*0.21);

c0 = 3*1E8;
VF = (2*pi.*f)./(c0 .* beta);

figure(2)
plot(f/1E6, VF, '-')
xlabel('Frequency (MHz)');
ylabel('Velocity Factor up/c0');
title('Velocity Factor');

xlim([0,300])
ylim([0.6,0.8])
yticks(0.6:0.05:0.8)
%% For Alpha
alpha = unwrap(real(gamma)) ./ (2*0.21); 

figure
plot(f/10^6, alpha, '-')
xlabel('Frequency (MHz)');
ylabel('Alpha (Np/m)');
title('Attenuation Constant');

xlim([0,1500])
ylim([0,0.5])
yticks(0:0.05:0.5)
