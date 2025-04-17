% Constants
c = 3e8; % Speed of light in vacuum (m/s)

% Antenna parameters
h = 0.25 / 2; % Half-length of the dipole in meters (125 mm)
a_values = [0.001/2, 4/2, 20/2] * 1e-3; % Radius values in meters (half diameters)

% Frequency range
freq = linspace(1e6, 2e9, 2000); % Frequencies from 1 MHz to 2 GHz

% Preallocate impedance arrays for different diameters
Zin = zeros(length(a_values), length(freq), 2); % Last dimension for R and X

% Calculate the input impedance for each diameter
for idx = 1:length(a_values)
    a = a_values(idx); % Current radius
    for f_idx = 1:length(freq)
        lambda = c / freq(f_idx); % Wavelength
        L = 2 * h; % Total length of the dipole
        
        % Radiation resistance and reactance formulas
        if L < lambda/2
            % Dipole shorter than half wavelength
            R = (pi^2 * (L/lambda)^2) / 6 * 120 * pi; % Approximate formula
            X = 60 * log((L/pi/a)) * (1 - (L/lambda)^2); % Approximate formula
        elseif L == lambda/2
            % Half-wave dipole
            R = 73; % Approximation for half-wave dipole in free space
            X = 0; % Reactance is zero at resonance
        else
            % Dipole longer than half wavelength, use empirical formulas or numerical methods
            R = 73 * (L / (lambda/2))^0.5; % Empirical scaling for longer dipoles
            X = 200 * (L/lambda - 0.5); % Rough empirical approximation
        end
        
        % Store the impedance values
        Zin(idx, f_idx, :) = [R, X];
    end
    
    % Plotting
    figure;
    subplot(1,2,1);
    plot(freq/1e6, squeeze(Zin(idx, :, 1)), 'b'); % Real part (R)
    hold on;
    plot(freq/1e6, squeeze(Zin(idx, :, 2)), 'r'); % Imaginary part (X)
    title(['Input Impedance for 2a = ', num2str(2*a*1e3), ' mm (Linear Scale)']);
    xlabel('Frequency (MHz)');
    ylabel('Impedance (\Omega)');
    legend('R', 'X');
    grid on;
    
    subplot(1,2,2);
    semilogy(freq/1e6, abs(squeeze(Zin(idx, :, 1))), 'b'); % Real part (R)
    hold on;
    semilogy(freq/1e6, abs(squeeze(Zin(idx, :, 2))), 'r--'); % Imaginary part (X)
    title(['Input Impedance for 2a = ', num2str(2*a*1e3), ' mm (Log Scale)']);
    xlabel('Frequency (MHz)');
    ylabel('Impedance (\Omega)');
    legend('R', 'X');
    grid on;
end
