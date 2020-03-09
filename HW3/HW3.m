% Maxwell Weil
% AMATH 482
% HW 3
clear; close all; clc;

% Setting number of camera and test number
numCams = 3;
dataSet = 1;

% Initializing name variables
camNames = strings(1,numCams);
vidNames = strings(1,numCams);

% Creating variable anmes
for i = 1:numCams
    camNames(i) = ['cam',num2str(i),'_',num2str(dataSet),'.mat'];
    vidNames(i) = ['vidFrames',num2str(i),'_',num2str(dataSet)];
end

% Initializing variables
numFrames = zeros(3,1);
camData = cell(1,3);
dataLen = zeros(1,3);
startIdx = zeros(1,3);

% Looping through each camera to find object location for each video
for i = 1:numCams
    
    % Loading in and setting current camera data
    load(camNames(i))
    currentVid = eval(vidNames(i));
    
    % Setting number of frames for current video
    numFrames(i) = size(currentVid, 4);
    
    % Initializing grayscale images
    grayImage = zeros(size(currentVid,1,2,4),'uint8');
    
    % Converting each frame of current video to grayscale
    for j = 1:numFrames(i)
        grayImage(:, :, j) = rgb2gray(currentVid(:, :, :, j));
    end
    
    % Selecting region of interest using rectangular area
    figure
    imshow(grayImage(:, :, 50))
    rectInfo = getrect;
    
    % Rounding pixel values to prevent any errors
    xCoor = round([rectInfo(1),rectInfo(1)+rectInfo(3)]);
    yCoor = round([rectInfo(2),rectInfo(2)+rectInfo(4)]);
    close
    
    % Initializing data and image difference matrices
    data = zeros(numFrames(i)-1,2);
    imageDiffs = zeros(size(currentVid,1,2),'uint8');
    
    % Looping through number of frames to find object location in each
    for k = 1:numFrames(i)-1
        
        % Calculating difference between images, in order to
        % highlights moving objects
        imageDiffs(:, :) = imabsdiff(grayImage(:, :, k), ...
            grayImage(:, :, k+1));
        
        % Binarizing image to reduce movement noise
        bw = imbinarize(imageDiffs(yCoor(1):yCoor(2), ...
            xCoor(1):xCoor(2)),0.3);

        % Find indices of nonzero values and averaging to determine
        % location of moving object
        [row, col] = find(bw);
        if ~isempty(row)
            data(k,2) = mean(row);
            data(k,1) = mean(col);
        end
        
        % Plotting binary image along with calculated object location
        imshow(bw)
        hold on
        plot(data(k,1),data(k,2),'rx','MarkerSize',20)
        hold off
        drawnow
    end
    close
    
    % Storing object location and rotating data from camera 3
    camData{i} = data;
    if i == 3
        camData{i}(:, [1 2]) = camData{i}(:, [2 1]);
    end
    
    dataLen(i) = length(camData{i});
    
    % Plotting all extracted object locations
    figure
    plot(camData{i}(:,1),camData{i}(:,2), 'r.', 'MarkerSize', 20)
    axis equal
    drawnow
    pause(3)
end

% Looping through camera to fill in any zero values and determine start
for i = 1:numCams
    
    % Setting data to look at for this loop iteration
    currentData = camData{i}(:,:);
    
    % Making any zero values NaN
    currentData(currentData==0) = NaN;
    
    % Finding number of existing data points
    nonNanValues = mean(sum(~isnan(currentData)));
    
    % Creating sample vector with length equal to data
    samples = linspace(1,nonNanValues,length(currentData));
    
    % Interpolating missing points at each NaN
    filledData = fillmissing(currentData,'linear','SamplePoints',samples);
    
    % Finding minimum in first handful of frames
    % And designating starting index for synchronization
    [val, idx] = min(filledData(10:40,2));
    startIdx(i) = idx;
    
    % Adjusting for new data length and saving filled data
    dataLen(i) = dataLen(i)-idx;
    camData{i} = filledData;
    
    
end

% Finding minimum data length to cut all data to
adjLen = min(dataLen);

% Initializing snapshot matrix
newData = zeros(numCams*2,adjLen-1);

% Looping through data and adjusting length
for i = 1:numCams
    
    % Selecting current data for this loop iteration
    currentData = camData{i}(:,:);
    
    % Filling in snapshot matrix with synchronized starts and lengths
    newData((2*i-1),:) = currentData(startIdx(i):startIdx(i)+adjLen-2,1);
    newData((2*i),:) = currentData(startIdx(i):startIdx(i)+adjLen-2,2);

end

% Setting mean of each row in snapshot matrix equal to zero
newData = newData-mean(newData,2);

% Calculating SVD
[U,S,V] = svd(newData,'econ');

% Finding sigma values from S matrix
sigmas = diag(S);

% Plotting sigma values and energy
figure
subplot(2,1,1)
plot(sigmas,'bo','Linewidth',2)
axis([0 length(sigmas) 0 max(sigmas)])
title('Sigma Values For Each Mode')
ylabel('\sigma Values')
xlabel('Mode')

subplot(2,1,2)
plot(sigmas.^2/sum(sigmas.^2),'bo','Linewidth',2)
axis([0 length(sigmas) 0 1])
title('Energy Accounted For in Each Mode')
ylabel('Energy Percent Contribution')
xlabel('Mode')
sgtitle(['Test ', num2str(dataSet), ' SVD Analysis'])
%%
% Calculating rank approximation for given rank number
rankNum = 1;
rankApprox = U(:,1:rankNum)*S(1:rankNum,1:rankNum)*V(:,1:rankNum)';

% Plotting data vs rank approximation for each camera in y directions
figure
for i = 1:size(newData,1)/2
    subplot(3,1,i)
    plot(newData(2*i,:), 'k', 'LineWidth', 3)
    hold on
    plot(rankApprox(2*i,:),'r.')
    hold off
    axis([0 length(newData) -100 100 ])
    title(['Camera ', num2str(i)])
    ylabel('Pixel Displacement (Y)', 'FontSize', 9)
    xlabel('Frame', 'FontSize', 9)
end
sgtitle(['Rank ', num2str(rankNum),...
    ' Approximations for Test ', num2str(dataSet)])
legend('Recorded Data', 'Approximated Data')

% Plotting data vs rank approximation for each camera in x directions
figure
for i = 1:size(newData,1)/2
    subplot(3,1,i)
    plot(newData(2*i-1,:), 'k', 'LineWidth', 3)
    hold on
    plot(rankApprox(2*i-1,:),'r.')
    hold off
    axis([0 length(newData) -100 100 ])
    title(['Camera ', num2str(i)])
    ylabel('Pixel Displacement (X)', 'FontSize', 9)
    xlabel('Frame', 'FontSize', 9)
end
sgtitle(['Rank ', num2str(rankNum),...
    ' Approximations for Test ', num2str(dataSet)])
legend('Recorded Data', 'Approximated Data')
