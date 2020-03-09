% Maxwell Weil
% AMATH 482
% HW 4

clear all; close all; clc;

% Requesting file numbers and names, assuming that song 
% Files are listed as ACDC1, ACDC2, Eminem1, Eminem2, etc.
num_files = input('How many files for each class (must be equal)?   ');
class1_name = input('First class file name? ', 's');
class2_name = input('Second class file name? ', 's');
class3_name = input('Third class file name? ', 's');

% Converting names to string arrays
class1_name = convertCharsToStrings(class1_name);
class2_name = convertCharsToStrings(class2_name);
class3_name = convertCharsToStrings(class3_name);

% Setting file names to read in
class_names = repelem([class1_name,class2_name,class3_name],num_files);
file_id = repmat(1:num_files,1,num_files);
file_names = strings(size(class_names));
for i=1:length(class_names)
    file_names(i) = strcat(class_names(i), num2str(file_id(i)));
end

% Setting number of training samples and test samples taken from each file,
% As well as duration of each sample in seconds
n_train = 20;
n_test = 10;
seconds = 5;

% Looping through files
for k = 1:length(file_names)
    
    % Reading in current file
    current_file = strcat(file_names(k),'.mp3');
    [y,fs] = audioread(current_file);

    % Downsampling audio and sampling rate
    y = downsample(y,4);
    fs = fs/4;

    % Calculating duration of song, and total possible samples that can be
    % Taken without any overlapping data
    duration = length(y)/fs;
    num_samples = floor(duration/seconds);
    
    % Randomly selecting nonoverlapping start points for samples
    rand_samples = randperm(num_samples-1,n_train+n_test);
    
    % Creating matrix of training data
    for i = 1:n_train
        train_data(:,i,k) = y(rand_samples(i)*seconds*fs: ...
            (rand_samples(i)+1)*seconds*fs-1);
    end
    
    % Creating matrix of test data
    for i = 1:n_test
        test_data(:,i,k) = y(rand_samples(i+n_train)*seconds*fs: ...
            (rand_samples(i+n_train)+1)*seconds*fs-1);
    end
end

% Reshaping data for spectrogram
train_data = reshape(train_data,size(train_data,1),[],1);
test_data = reshape(test_data,size(test_data,1),[],1);

% Looping through training data and computing spectrogram
for j = 1:size(train_data,2)
    spec = spectrogram(train_data(:,j)',gausswin(1000), ...
        560, 246, fs, 'yaxis');
    train_spec(:,j) = reshape(abs(spec),[],1);
end

% Computing spectrogram for test data, parameters were chosen for optimal
% Performance given the type and size of data
for j = 1:size(test_data,2)
    spec = spectrogram(test_data(:,j)',gausswin(1000), ...
        560, 246, fs, 'yaxis');
    test_spec(:,j) = reshape(abs(spec),[],1);
end

% Setting number of features to look at in training data
feature = 20;

% Computing SVD of training data, making feature space
[U,S,V] = svd(train_spec, 'econ');
music = S*V';
U = U(:,1:feature);

% Finding number of training data for each class (should be equal)
n_class1 = num_files*n_train;
n_class2 = num_files*n_train;
n_class3 = num_files*n_train;

% Finding columns of music matrix pertaining to each class
class1 = music(1:feature,1:n_class1);
class2 = music(1:feature,n_class1+1:n_class1+n_class2);
class3 = music(1:feature,n_class1+n_class2+1:n_class1+n_class2+n_class3);

% Finding average values of PCA projection for each class, and overall
avg_class1 = mean(class1,2);
avg_class2 = mean(class2,2);
avg_class3 = mean(class3,2);
avg_tot = mean(music(1:feature,:),2);

% Computing in class variances
Sw = 0;
for k=1:n_class1
    Sw = Sw + (class1(:,k)-avg_class1)*(class1(:,k)-avg_class1)';
end

for k=1:n_class2
    Sw = Sw + (class2(:,k)-avg_class2)*(class2(:,k)-avg_class2)';
end

for k=1:n_class2
    Sw = Sw + (class3(:,k)-avg_class3)*(class3(:,k)-avg_class3)';
end

% Computing between class variances
Sb_class1 = (avg_class1-avg_tot)*(avg_class1-avg_tot)';
Sb_class2 = (avg_class2-avg_tot)*(avg_class2-avg_tot)';
Sb_class3 = (avg_class3-avg_tot)*(avg_class3-avg_tot)';
Sb = Sb_class3 + Sb_class2 +Sb_class1;

% Computing LDA
[V2,D] = eig(Sb,Sw);

% Finding the two nonzero eigenvalues and their corresponding vectors
[~, ind] = sort(abs(diag(D)), 'descend');
w = V2(:,ind(1:2));
w = w/norm(w,2);

% Projecting each class onto both eigenvectors
proj_class1 = w'*class1; 
proj_class2 = w'*class2;
proj_class3 = w'*class3;

% Finding the center of each class cluster
centroids = [mean(proj_class1(1,:)), mean(proj_class1(2,:));
    mean(proj_class2(1,:)), mean(proj_class2(2,:));
    mean(proj_class3(1,:)), mean(proj_class3(2,:))];

% Establishing number of test samples for each class
n_test_samp = n_test*num_files;

% Projecting test data onto trained PCA and LDA found previously
test_proj_PCA = U'*test_spec;
test_proj_LDA = w'*test_proj_PCA;

% Initializing prediction vector
predictions = zeros(1,n_test*length(class_names));

% Looping through LDA projections for test sample
for i = 1:length(test_proj_LDA)
    
    % Initializing distance vector
    dist = zeros(1,3);
    
    % Looping through each class's center
    for j = 1:size(centroids,1)
        
        % Finding distance of LDA projection from each class center
        dist(j) = pdist2(test_proj_LDA(:,i)',centroids(j,:));
    end
    
    % Computing the closest class center and saving as prediction
    [M,I] = min(dist);
    predictions(i) = I;
end

% Creating vector of true classes for each test sample
truths = repelem([1:3],n_test_samp);

% Finding correct and incorrect predictions
correct = predictions==truths;
incorrect = ~correct;

% Computing percent correct
percent_correct = sum(correct)/length(correct)*100;
%%
% Creating plotting vectors for correct prediction in each class 
% Each row is 1 for all correct predictions for that class,
% And 0 for incorrect predictions or other classes
correct_plot = logical([correct(1:n_test_samp) zeros(1,2*n_test_samp);
    zeros(1,n_test_samp) correct(n_test_samp+1:2*n_test_samp) ...
        zeros(1,n_test_samp);
    zeros(1,2*n_test_samp) correct(2*n_test_samp+1:3*n_test_samp)]);

% Plotting training and test data
% As well aslong with correct and incorrect predictions
plot(proj_class1(1,:), proj_class1(2,:), '.m')
hold on
plot(proj_class2(1,:), proj_class2(2,:), '.b')
plot(proj_class3(1,:), proj_class3(2,:), '.k')
plot(test_proj_LDA(1,correct_plot(1,:)), ... 
    test_proj_LDA(2,correct_plot(1,:)), 'mo')
plot(test_proj_LDA(1,correct_plot(2,:)), ...
    test_proj_LDA(2,correct_plot(2,:)), 'bo')
plot(test_proj_LDA(1,correct_plot(3,:)), ...
    test_proj_LDA(2,correct_plot(3,:)), 'ko')
plot(test_proj_LDA(1,incorrect),test_proj_LDA(2,incorrect), 'rx')

xlabel('LDA Projection on First Eigenvector')
ylabel('LDA Projection on Second Eigenvector')
title(["Projections of Training and Test Data on 2D LDA Space", ...
    strcat("Prediction Accuracy of ", num2str(percent_correct), "%")])
legend(strcat(class1_name," Training Data"), ...
    strcat(class2_name," Training Data"), ...
    strcat(class3_name," Training Data"), ...
    strcat(class1_name," Test Data (Correct)"), ...
    strcat(class2_name," Test Data (Correct)"), ...
    strcat(class3_name," Test Data (Correct)"), ...
    "Incorrect",'Location','best')

% Computing and printing number of misclassified points from each class
class1_in = sum(~correct_plot(1,:))-2*n_test_samp;
class2_in = sum(~correct_plot(2,:))-2*n_test_samp;
class3_in = sum(~correct_plot(3,:))-2*n_test_samp;
disp(strcat(class1_name, " misclassified = ", num2str(class1_in)))
disp(strcat(class2_name, " misclassified = ", num2str(class2_in)))
disp(strcat(class3_name, " misclassified = ", num2str(class3_in)))

%% Producing sample plots and spectrograms

% Creating time vector for time series plotting
time = 1/fs:1/fs:seconds;

% Choosing last sample from each class
class1_samp = train_data(:,n_class1);
class2_samp = train_data(:,n_class1+n_class2);
class3_samp = train_data(:,n_class1+n_class2+n_class3);

% Plotting time series of sample from each class
subplot(3,1,1)
plot(time, class1_samp, 'm')
title(strcat(class1_name, " Time Series Sample"))
axis([0 5 -1 1])
xlabel('Time')
ylabel('Amplitude')
subplot(3,1,2)
plot(time, class2_samp, 'b')
title(strcat(class2_name, " Time Series Sample"))
axis([0 5 -1 1])
xlabel('Time')
ylabel('Amplitude')
subplot(3,1,3)
plot(time, class3_samp, 'k')
title(strcat(class3_name, " Time Series Sample"))
axis([0 5 -1 1])
xlabel('Time')
ylabel('Amplitude')

% Plotting spectrogram sample from each class
figure
subplot(3,1,1)
spectrogram(class1_samp',gausswin(1000), 560, 246, fs, 'yaxis');
title(strcat(class1_name, " Spectrogram Sample"))
subplot(3,1,2)
spectrogram(class2_samp',gausswin(1000), 560, 246, fs, 'yaxis');
title(strcat(class2_name, " Spectrogram Sample"))
subplot(3,1,3)
spectrogram(class3_samp',gausswin(1000), 560, 246, fs, 'yaxis');
title(strcat(class3_name, " Spectrogram Sample"))
colormap bone