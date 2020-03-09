% Maxwell Weil
% AMATH 482
% HW2
clear; close all; clc;

% Loading sample music and sampling rate, reorienting vector
load handel
music = y';

% Setting sample number, duration of music, time vector,
% And unshifted/shifted frequency vectors 
samples = length(music);
duration = samples/Fs;
time = [(1:samples)/Fs];
freqs = (1/duration)*[0:(samples-1)/2 -(samples-1)/2:-1];
shifted_freqs=fftshift(freqs);


% Defining constant for plotting static windows
a_ex = 10;
width_ex = 3;
sigma_ex = 10;
tau = 4;

% Defining equations for gabor, shannon, and mexican hat windows
gabor_window=exp(-a_ex*(time-tau).^2);
shannon_window = heaviside(time-(tau-width_ex/2))...
    -heaviside(time-(tau+width_ex/2));
mexican_hat_window = (1-sigma_ex*(time-tau).^2)...
    .*exp(-1*sigma_ex*(time-tau).^2/2);

% Creating plots of static windows and filtered signal in
% Frequency and time domains
figure

% Gabor window subplots
subplot(3,3,1)
plot(time, gabor_window, 'r', 'LineWidth', 2)
axis([0 duration -0.5 1.5])
title('Gabor Window', 'FontSize', 6)
xlabel('Time (s)', 'FontSize', 6)
ylabel('Amplitude', 'FontSize', 6)
subplot(3,3,2)
plot(time, gabor_window, 'r', 'LineWidth', 2)
hold on
plot(time, gabor_window.*music, 'k')
axis([0 duration -1 1.5])
title('Gabor Filtered Signal', 'FontSize', 6)
xlabel('Time (s)', 'FontSize', 6)
ylabel('Amplitude', 'FontSize', 6)
subplot(3,3,3)
plot(shifted_freqs,abs(fftshift(fft(gabor_window.*music)))...
    /max(abs(fft(gabor_window.*music))),'m');
axis([-Fs/2 Fs/2 0 1])
title('Gabor Filtered Signal (Frequency)', 'FontSize', 6)
xlabel('Frequency (Hz)', 'FontSize', 6)
ylabel('Relative Amplitude', 'FontSize', 6)

% Shannon window subplots
subplot(3,3,4)
plot(time, shannon_window, 'r', 'LineWidth', 2)
axis([0 duration -0.5 1.5])
title('Shannon Window', 'FontSize', 6)
xlabel('Time (s)', 'FontSize', 6)
ylabel('Amplitude', 'FontSize', 6)
subplot(3,3,5)
plot(time, shannon_window, 'r', 'LineWidth', 2)
hold on
plot(time, shannon_window.*music, 'k')
axis([0 duration -1 1.5])
title('Shannon Filtered Signal', 'FontSize', 6)
xlabel('Time (s)', 'FontSize', 6)
ylabel('Amplitude', 'FontSize', 6)
subplot(3,3,6)
plot(shifted_freqs,abs(fftshift(fft(shannon_window.*music)))...
    /max(abs(fft(shannon_window.*music))),'m');
axis([-Fs/2 Fs/2 0 1])
title('Shannon Filtered Signal (Frequency)', 'FontSize', 6)
xlabel('Frequency (Hz)', 'FontSize', 6)
ylabel('Relative Amplitude', 'FontSize', 6)

% Mexican hat window subplots
subplot(3,3,7)
plot(time, mexican_hat_window, 'r', 'LineWidth', 2)
axis([0 duration -0.5 1.5])
title('Mexican Hat Window', 'FontSize', 6)
xlabel('Time (s)', 'FontSize', 6)
ylabel('Amplitude', 'FontSize', 6)
subplot(3,3,8)
plot(time, mexican_hat_window, 'r', 'LineWidth', 2)
hold on
plot(time, mexican_hat_window.*music, 'k')
axis([0 duration -1 1.5])
title('Mexican Hat Filtered Signal', 'FontSize', 6)
xlabel('Time (s)', 'FontSize', 6)
ylabel('Amplitude', 'FontSize', 6)
subplot(3,3,9)
plot(shifted_freqs,abs(fftshift(fft(mexican_hat_window.*music)))...
    /max(abs(fft(mexican_hat_window.*music))),'m');
axis([-Fs/2 Fs/2 0 1])
title('Mexican Hat Filtered Signal (Frequency)', 'FontSize', 6)
xlabel('Frequency (Hz)', 'FontSize', 6)
ylabel('Relative Amplitude', 'FontSize', 6)
%% Effect of Different Window Sizes on Spectrogram

% Defining window size constants for sliding gabor window
a_list = [1, 10, 100, 10000];

% Creating time steps to move window along
tslide=0:0.01:duration;

% Initializing spectrogram matrix 
gab_spec_music = zeros(length(tslide),samples);

% Looping through different sized windows for gabor-filtered spectrogram
figure
for i = 1:length(a_list)
    
    % Looping through time step vector
    for j=1:length(tslide)
        
        % Creating gabor window at time step
        gabor_window=exp(-a_list(i)*(time-tslide(j)).^2);
        
        % Filtering signal in time
        gab_filt_music=gabor_window.*music; 
        
        % Transforming signal to frequency
        freq_gab_filt_music=fft(gab_filt_music); 
        
        % Adding frequency data to spectrogram at each time step
        gab_spec_music(j,:) = fftshift(abs(freq_gab_filt_music));
        
    end
        % Plotting spectrograms for each sized window
        subplot(2,2,i)
        pcolor(tslide,shifted_freqs,gab_spec_music.')
        ylim([0 max(freqs)])
        title(['Gabor Filtered Spectrogram with a = ', num2str(a_list(i))])
        ylabel('Frequency (Hz)')
        xlabel('Time (s)')
        shading interp 
        colormap hot
end
%% Effect of Step Size on Spectrogram

% Defining window size constant
a = 100;

% Creating varying time steps to move window along
step_list = [0.05, 0.1, 0.5, 1];

% Looping through different sized time steps for gabor-filtered spectrogram
figure
for i = 1:length(step_list)
    
    % Setting new time step size
    tslide=0:step_list(i):duration;
    
    % Initializing spectrogram matrix 
    gab_spec_music = zeros(length(tslide),samples);
    
    % Looping through time step vector
    for j=1:length(tslide)
        
        % Creating gabor window at time step
        gabor_window=exp(-a*(time-tslide(j)).^2);
        
        % Filtering signal in time
        gab_filt_music=gabor_window.*music; 
        
        % Transforming signal to frequency
        freq_gab_filt_music=fft(gab_filt_music); 
        
        % Adding frequency data to spectrogram at each time step
        gab_spec_music(j,:) = fftshift(abs(freq_gab_filt_music));
        
    end
        % Plotting spectrograms for each sized window
        subplot(2,2,i)
        pcolor(tslide,shifted_freqs,gab_spec_music.')
        ylim([0 max(freqs)])
        title({'Gabor Filtered Spectrogram'; 'with time step = '; ...
              num2str(step_list(i)))
        ylabel('Frequency (Hz)')
        xlabel('Time (s)')
        shading interp 
        colormap hot
end
%% Effect of Different Window Types on Spectrogram

% Defining optimized constants for shannon and mexican hat windows
a = 100;
width = 0.1;
sigma = 100;


% Initializing spectrogram matrices for each window type
gab_spec_music = zeros(length(tslide),samples);
shan_spec_music = zeros(length(tslide),samples);
mex_spec_music = zeros(length(tslide),samples);

% Looping through time steps, with optimized
% Window sizes, for spectrogram comparison
figure
for j=1:length(tslide)
    
   % Creating each window at specified time step
   gabor_window=exp(-a*(time-tslide(j)).^2);
   shannon_window = heaviside(time-(tslide(j)-width/2))-heaviside(time-(tslide(j)+width/2));
   mexican_hat_window = (1-sigma*(time-tslide(j)).^2).*exp(-1*sigma*(time-tslide(j)).^2/2);
   
   % Filtering and transforming signal with gabor window
   gab_filt_music=gabor_window.*music; 
   freq_gab_filt_music=fft(gab_filt_music);
   
   % Filtering and transforming signal with shannon window
   shan_filt_music=shannon_window.*music; 
   freq_shan_filt_music=fft(shan_filt_music);
   
   % Filtering and transforming signal with mexican hat window
   mex_filt_music=mexican_hat_window.*music; 
   freq_mex_filt_music=fft(mex_filt_music);
   
   % Adding frequency data to each spectrogram at each time step
   gab_spec_music(j,:) = fftshift(abs(freq_gab_filt_music));
   shan_spec_music(j,:) = fftshift(abs(freq_shan_filt_music));
   mex_spec_music(j,:) = fftshift(abs(freq_shan_filt_music));
end

% Plotting spectrograms for each window type
% Gabor spectrogram
subplot(1,3,1)
pcolor(tslide,shifted_freqs,gab_spec_music.')
ylim([0 max(freqs)])
title({'Gabor Filtered'; 'Spectrogram'})
ylabel('Frequency (Hz)')
xlabel('Time (s)')
shading interp 
colormap hot

% Shannon spectrogram
subplot(1,3,2)
pcolor(tslide,shifted_freqs,shan_spec_music.')
ylim([0 max(freqs)])
title({'Shannon Filtered'; 'Spectrogram'})
ylabel('Frequency (Hz)')
xlabel('Time (s)')
shading interp 
colormap hot

% Mexican hat spectrogram
subplot(1,3,3)
pcolor(tslide,shifted_freqs,mex_spec_music.')
ylim([0 max(freqs)])
title({'Mexican Hat Filtered'; 'Spectrogram'})
ylabel('Frequency (Hz)')
xlabel('Time (s)')
shading interp 
colormap hot

%% PROBLEM 2, Analyzing A Music Score
clear; close all; clc;

% Reading in audio data and sampling frequency
[y,Fs] = audioread('music2.wav');

% Downsampling audio for fast computing, note that Fs is also adjusted
y = downsample(y,4);
Fs = Fs/4;

% Reorienting audio vector, setting sample number, duration of music,
% Time vector, and unshifted/shifted frequency vectors 
song = y';
samples = length(song);
duration = samples/Fs;
time = [(1:samples)/Fs];
freqs = (1/duration)*[0:samples/2-1 -samples/2:-1];
shifted_freqs=fftshift(freqs);

% Setting selected filter constants and time step vector
a = 100;
a2 = 0.01;
tslide=0:0.1:duration;

% Initializing spectrogram, maximum amplitude, and tone
% At maximum amplitude matrices
spec_song = zeros(length(tslide),samples);
max_amplitude = zeros(length(tslide),1);
tone_at_max = zeros(length(tslide),1);

% Looping through time steps
for j=1:length(tslide)
    
    % Creating window and filtering audio
    gabor_filter=exp(-a*(time-tslide(j)).^2);
    filt_song=gabor_filter.*song; 
    
    % Transforming to frequency domain
    filt_freq_song=fft(filt_song);
    
    % Finding the maximum amplitude, and the frequency at which it occurs
    [high, idx] = max(abs(filt_freq_song));
    tone_at_max(j) = freqs(idx);
    max_amplitude(j) = high;
    
    % Creating filter in frequency domain around center frequency
    overtone_filter = exp(-a2*(freqs-freqs(idx)).^2);
    
    % Filtering out overtones around center frequency
    pure_freq_song = overtone_filter.*filt_freq_song;
    
    % Adding frequency data to each spectrogram at each time step
    spec_song(j,:) = fftshift(abs(filt_freq_song));
end

% Plotting spectrogram
figure
pcolor(tslide,shifted_freqs,spec_song.')
ylim([0 2400])
title('Recorder Spectrogram')
ylabel('Frequency (Hz)')
xlabel('Time (s)')
shading interp
colormap hot

% Finding where the maximum amplitudes are over a specified threshold,
% In order to determine true notes are being played
threshold_logic = max_amplitude>mean(max_amplitude)/2;

% Removing repeated values, by finding difference between adjacent elements
high_low_shift = diff(threshold_logic) < 0;

% Locating corresponding frequencies for each note
notes_freq = tone_at_max(high_low_shift);

% Initializing text string for notes
notes_text = strings(length(notes_freq),1);

% Looping through frequncies and determining alphabetical note,
% This process has been shortened to only search for a few specific notes
for i = 1:length(notes_freq)
    if notes_freq(i) < 1100 && notes_freq(i) > 980
        notes_text(i) = 'B';
    elseif notes_freq(i) < 980 && notes_freq(i) > 850
        notes_text(i) = 'A';
    elseif notes_freq(i) < 850 && notes_freq(i) > 750
        notes_text(i) = 'G';
    elseif notes_freq(i) < 330 && notes_freq(i) > 300
        notes_text(i) = 'E';
    elseif notes_freq(i) < 300 && notes_freq(i) > 280
        notes_text(i) = 'D';
    else
        notes_text(i) = 'C';
    end
end

% Creating vector to plot notes
note_scaling = linspace(1, 100*length(notes_freq), length(notes_freq));

% Plotting frequencies and corresponding alphabetical note, keep in mind
% That the x-axis is arbitrary, and notes are not actually equally spaced
figure
plot(note_scaling,notes_freq, 'ko', 'MarkerSize', 10)
text(note_scaling,notes_freq+20, notes_text)
axis([(min(note_scaling)-100) (max(note_scaling)+100)...
    (min(notes_freq)-300) (max(notes_freq)+300)])
title('Reproduced Music Score for Piano')
xlabel('Relative Time (a.u.)')
ylabel('Frequency (Hz)')