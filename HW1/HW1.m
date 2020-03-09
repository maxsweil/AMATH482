% Maxwell Weil
% AMATH 482
% 1/24/2020

clear; close all; clc;

% Read in data
load Testdata
L=15; % Spatial domain
n=64; % Fourier modes
trials = 20; % Total number of trials over time

% Setting dimensions for spatial domain
x2=linspace(-L,L,n+1);
x=x2(1:n);
y=x;
z=x;

% Setting frequencies for frequency domain
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1];

% Creating full grid of spatial and frequency axes
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(k,k,k);

% Initializing average frequency matrix [size of data]
avg_freq_Un = zeros(n,n,n);

% Looping across all time and summing frequencies with reshaped data
for j=1:trials
    Un(:,:,:)= reshape(Undata(j,:),n,n,n);
    avg_freq_Un = avg_freq_Un + fftn(Un);
end

% Calculating average frequency of signal over all time
avg_freq_Un = abs(avg_freq_Un)/trials;

% Finding maximum frequency component (center frequency) and index
[maxval, idx] = max(avg_freq_Un,[],'all','linear');

% Converting linear index to 64x64x64 indices
[row, col, lay] = ind2sub([n n n], idx);

% Finding frequency center in 3D frequency domain
Kx0 = Kx(row, col, lay);
Ky0 = Ky(row, col, lay);
Kz0 = Kz(row, col, lay);

% Setting tau constants for spectral filter
taux = 0.2;
tauy = 0.1;
tauz = 0.3;

% Creating gaussian spectral filter
filter = exp(-taux*(Kx-Kx0).^2-tauy*(Ky-Ky0).^2-tauz*(Kz-Kz0).^2);

% Plot filter in frequency domain
figure(1)
isosurface(fftshift(Kx),fftshift(Ky),fftshift(Kz),fftshift(filter), 0.9)
title('Filter in Frequency Domain', 'FontSize', 18)
axis([1 3 -3 1 -1 1])
xlabel('Kx Direction'); ylabel('Ky Direction'); zlabel('Kz Direction')

% Initializing marble locations matrix [dimensions x time points]
marb_loc = zeros(3,trials);

figure(2)
for j = 1:trials
    
    % Reshaping data at each time point
    Un(:,:,:)=reshape(Undata(j,:),n,n,n);
    
    % Change data into frequency domain
    avg_freq_Un = fftn(Un);
    
    % Applying spectral filter to frequency data
    filt_freq_Un = filter .* avg_freq_Un; 
    
    % Transforming data back into spacial domain
    filt_space_Un = ifftn(filt_freq_Un);
    
    % Finding maximum spatial point, converting to true axes
    [maxval2, idx2] = max(filt_space_Un,[],'all','linear');
    [row, col, lay] = ind2sub([n n n], idx2);
    marb_loc(1,j) = X(row, col, lay);
    marb_loc(2,j) = Y(row, col, lay);
    marb_loc(3,j) = Z(row, col, lay);
    
    % Plotting spatial contour at each time, animated
    isosurface(X,Y,Z,abs(filt_space_Un),0.4)
    title('Filtered Spatial Ultrasound Over Time', 'FontSize', 18)
    axis([-15 15 -15 15 -15 15])
    xlabel('X Direction'); ylabel('Y Direction'); zlabel('Z Direction')
    hold on; drawnow; pause(0.2);
end


% Plotting center of marble at each time, animated
figure(3)
for j = 1:trials
    plot3(marb_loc(1,1:j),marb_loc(2,1:j),marb_loc(3,1:j), ...
        'bo-','MarkerSize',10,'MarkerFaceColor','c')
    title('Marble Center Location Path', 'FontSize', 18)
    axis([-15 15 -15 15 -15 15])
    xlabel('X Direction'); ylabel('Y Direction'); zlabel('Z Direction')
    hold on; drawnow; pause(0.2);
end

% Display coordinates of final position
disp(['Send acoustic pulse to X = ', num2str(marb_loc(1,20)), ...
    ' Y = ', num2str(marb_loc(2,20)), ' Z = ', num2str(marb_loc(3,20))])