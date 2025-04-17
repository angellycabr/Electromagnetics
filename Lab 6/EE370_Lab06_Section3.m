clear;
clc;
close all;
%% Section 3
radii = [0.001, 4, 20] / 1000; % Conversion from mm to m
h = 250 / 2000; % Height in meters

% Specify z range
npoints = 500;
zstrt = 0;
zstop = h;
zstep = (zstop - zstrt) / (npoints - 1);
z = -zstop:zstep:zstop;

% Specify frequencies
nfreq = 400;
fstrt = 0 * 10^9;
fstop = 2 * 10^9;
fstep = (fstop - fstrt) / (nfreq - 1);
freq = fstrt:fstep:fstop;
lambda = (3 * 1E8) ./ freq;

% Define parameters
eta = 120 * pi;
Im = 1;

% Loop over each radius
for r_idx = 1:length(radii)
    r = radii(r_idx);
    
    R0 = sqrt(r^2 + z.^2);
    R1 = sqrt(r^2 + (z - h).^2);
    R2 = sqrt(r^2 + (z + h).^2);
    
    Zi = zeros(1, length(freq));
    
    for i = 1:length(freq)
        beta = (2 * pi) ./ lambda(i);
        
        E_z = ((1j * eta) / (4 * pi)) .* (((2 * exp(-1j .* beta * R0) * cos(beta * h)) ./ R0) - (exp(-1j .* beta * R1) ./ R1) - (exp(-1j .* beta * R2) ./ R2)); 
        H = (-1j * Im) / (4 * pi) * (1 / r) .* (2 * exp(-1j .* beta * R0) .* cos(beta * h) - exp(-1j .* beta * R1) - exp(-1j .* beta * R2));
        
        %In = E_z .* conj(H);
        %Zi(i) = ((-2 * pi * r) / (sin(beta * h)).^2) * trapz(z, In); %Eq. 25

        In = E_z .* sin(beta .* (h - abs(z)));
        Zi(i)=(-1 / (sin(beta .* h)).^2) .* trapz(z,In); %Eq. 30
    end

    % Plotting
    figure;
    title(['Dipoles Input Impedance, Zi=R+jX (2a=' num2str(radii(r_idx)*1000) ' mm)']);
    yyaxis left
    plot(freq / 1E6, real(Zi));
    xlim([0, 2000])
    xticks(0:100:2000)
    ylim([0, 400])
    yticks(0:100:400)
    xlabel('Frequency (MHz)')
    ylabel('R (in ohms)');
    
    hold on
    yyaxis right
    plot(freq / 1E6, imag(Zi));
    ylim([-1000, 1000]);
    yticks(-1000:500:1000)
    ylabel('X (in ohms)');

    % Natural Log
    figure;
    title(['Dipoles Input Impedance, Zi=R+jX (2a=' num2str(radii(r_idx)*1000) ' mm)']);
    yyaxis left
    plot(freq / 1E6, log(real(Zi)));
    xlim([0, 2000])
    xticks(0:500:2000)
    ylim([-4, 8])
    yticks(-4:2:8)
    xlabel('Frequency (MHz)')
    ylabel('log10 R (in ohms)');
    
    hold on
    yyaxis right
    plot(freq / 1E6, log(abs(imag(Zi))));
    ylim([0, 8]);
    yticks(0:2:8)
    ylabel('log10 X (in ohms)');
end
%% 3.6
r = 0.4/1000;
h = 0.120; %250 / 2000;

% Specify z range
npoints = 700;
zstrt = 0;
zstop = h;
zstep = (zstop - zstrt) / (npoints - 1);
z = -zstop:zstep:zstop;

% Specify frequencies
nfreq = 600;
fstrt = 0 * 10^9;
fstop = 2 * 10^9;
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
    
for i = 1:length(freq)
    beta = (2 * pi) ./ lambda(i);
    E_z = ((1j * eta) / (4 * pi)) .* (((2 * exp(-1j .* beta * R0) * cos(beta * h)) ./ R0) - (exp(-1j .* beta * R1) ./ R1) - (exp(-1j .* beta * R2) ./ R2)); 
    H = (-1j * Im) / (4 * pi) * (1 / r) .* (2 * exp(-1j .* beta * R0) .* cos(beta * h) - exp(-1j .* beta * R1) - exp(-1j .* beta * R2));
        
    %In = E_z .* conj(H);
    %Zi(i) = ((-2 * pi * r) / (sin(beta * h)).^2) * trapz(z, In); %Eq. 25

    In = E_z .* sin(beta .* (h - abs(z)));
    Zi(i)=(-1 / (sin(beta .* h)).^2) .* trapz(z,In); %Eq. 30

end

figure;
title(['Dipoles Input Impedance, Zi=R+jX (2h=' num2str(2*h*1000) ' mm)']);
yyaxis left
plot(freq / 1E6, 20*log10(real(Zi)));
xlim([550, 650])
xticks(550:50:650)
ylim([-100, 100])
yticks(-100:50:100)
xlabel('Frequency (MHz)')
ylabel('R (in ohms)');
    
hold on
yyaxis right
plot(freq / 1E6, 20*log10(abs(imag(Zi))));
ylim([-100, 100])
yticks(-100:50:100)
ylabel('X (in ohms)');