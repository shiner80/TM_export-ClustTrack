%% Matlab Script for analysis of TRACKMATE export XLM files, quantification of single-cluster MSD and measurement of cluster size

% This routine reads in xlm files exported from the TRACKMATE Fiji plug-in
% and the corresponding time-lapse movie, and calculates MSDs and particle
% sizes for each of the tracked clusters.
%____________________________________________
%
% 2014. Davide Mazza. San Raffaele Scientific Institute. Milan. Italy
%____________________________________________



%% Prepare workspace [Clears the entire content of the workspace]

clear; clc;

%% SET HERE the parameters for the analysis:


% TAG Measure Size of the spots
% If this Tag is set to 1, the routines go back to the images to calculate
% the size of the spots


MeasureSize = 1;



%Number of channels Set if the movie has one or two channels
%Two-color images must be saved as RGB stacks.
nChannels = 1;


pixelSize = 0.18674;     % Size of the pixel in µm
tStep = 0.5;               % Time between two frames of the movie in sec
nMSDdisplay = 20;       % Number of MSD time points to be calculated and displayed
nMSDfit = 20;            % Number of MSD time points to be fit with anomalous diffusion

fitWindow = 10;          % Size of window for fitting gaussians to the spots; [pixels]



% Calculate timepoints for MSD calculation
AA_tlist = [tStep: tStep: tStep*nMSDdisplay];



%% Read-In Track Data
[FileNameData,PathNameData]=uigetfile('.xml', 'Read In TrackMate File');
[tracks, md] = importTrackMateTracks([PathNameData,FileNameData]);


%% Read-In Tif File
[FileNameImg,PathNameImg]=uigetfile('.tif', 'Read In RGB Stack');
Img = TIFread([PathNameImg,FileNameImg]);
[height, width] = size(Img(1).data(:,:,1));  % Measure size of image


%% Select a region containing the cell for the calculation of the centroid

% The routine will display in RED the first image of the time-series and in
% GREEN the last image of the series

figure;
image(cat(3,Img(1).data(:,:,1),Img(end).data(:,:,1),Img(end).data(:,:,3)));
rect = round(getrect);
title('First and last frame of the movie - Select a region containing the entire cell')
Mask = zeros(height,width);
Mask(rect(2): rect(2) + rect(4), rect(1): rect(1)+rect(3)) = 1;
close;


% Calculate Centroid for each of the frames
for i = 1:length(Img)
    ImgMask = double(Img(i).data(:,:,1)).*Mask;
    
    X_hist=sum(ImgMask,1);
    Y_hist=sum(ImgMask,2);
    X=1:width; Y=1:height;
    AA_cent(i,1) =sum(X.*X_hist)/sum(X_hist);
    AA_cent(i,2) = sum(Y'.*Y_hist)/sum(Y_hist);
    
end

% Convert the centroid to physical units
AA_cent(:,1) = (AA_cent(:,1) - AA_cent(1,1))*pixelSize;
AA_cent(:,2) = (AA_cent(:,2) - AA_cent(1,2))*pixelSize;

%% Process Tracks
AA_OUT_MSD = [];
AA_OUT_MSDerr = [];
AA_OUT_COEF = [];

for i = 1:length(tracks)
    track = tracks{i};
    track = sortrows(tracks{i},1);     %sort rows
    
    %find and fill gaps in tracks
    idx_gaps = find(track(2:end,1)- track(1:end-1,1) > 1); % find gaps
    if ~isempty(idx_gaps)
        
        for j = 1:length(idx_gaps);
            k = idx_gaps(j);
            N_steps = track(k+1,1)-track(k,1) + 1;
            trackTemp = [track(k,1):track(k+1,1)]';
            trackTemp(:,2) = linspace(track(k,2),track(k+1,2), N_steps);
            trackTemp(:,3) = linspace(track(k,3),track(k+1,3), N_steps);
            trackTemp(:,4) = 0;
            trackTemp = trackTemp(2:end-1,:);
            track = [track;trackTemp];
        end
        track = sortrows(track);
    end
    
    tracks{i}= track;
    
    % CORRECT TRACKS FOR THE CELL MOVEMENT BY SUBTRACTING THE CENTROID
    % POSITION TO THE TRACJS.
    
    trackLim = min(track(:,1)) : max(track(:,1));
    track(:,2:3) = track(:,2:3) - AA_cent(trackLim + 1,:);
    
    tracks_corr{i} = track;
    
    
    
    %% Calculate mean squared displacements for Tracks and Fit them with anomalous diffusion fit
    
    
    testTrack = [];
    
    testTrack(:,1) = tracks_corr{i}(:,2);
    testTrack(:,2) = tracks_corr{i}(:,3);
    testTrack(:,3) = tracks_corr{i}(:,1);
    testTrack(:,4) = 1;
    
    [MSD, MSD_COEF] = calculateMSD_AD(testTrack, AA_tlist, nMSDfit, 0);
    
    AA_OUT_MSD(:,i) = MSD(:,2);
    AA_OUT_MSDerr(:,i) = MSD(:,3);
    AA_OUT_COEF(1:2,i) = MSD_COEF;
    
    
    
    
    
    %% Measure Size and intensity of each of the particles in both channels
    if MeasureSize
        frameId = tracks{i}(1,1) + 1;
        xcoord = tracks{i}(1,2)/pixelSize +1;         % initial coordinate of the particle
        ycoord = tracks{i}(1,3)/pixelSize +1;
        
        ImRed =   Img(frameId).data(:,:,1);
        
        
        % for each spot fit gaussian with fixed center
        
        % Cut sub image
        
        xcoord_round = round(xcoord);
        ycoord_round = round(ycoord);
        
        ysub = (ycoord_round-fitWindow/2:ycoord_round+fitWindow/2);
        xsub = (xcoord_round-fitWindow/2:xcoord_round+fitWindow/2);
        
        
        ImRedsub = double(ImRed(ysub,xsub));
        
        coordinates = {xsub, ysub};
        center = [xcoord,ycoord];
        
        
        % Ft Gaussian
        [parFitRed, ssr] = gauss_fit2D_cov_fixedC(ImRedsub,coordinates, center);
        FitImgRed = gaussFun2D_cov_fixedC(parFitRed, coordinates);
        
        
        
        % Export measurements on particle size and total intensity
        AA_OUT_COEF(3,i) = parFitRed(2);
        AA_OUT_COEF(4,i) = (parFitRed(1)*parFitRed(2)*sqrt(2*pi))^2;
        
        
        if nChannels == 2;
            ImGreen = Img(frameId).data(:,:,2);
            ImGreensub = double(ImGreen(ysub,xsub));
            [parFitGreen, ssr] = gauss_fit2D_cov_fixedC(ImGreensub,coordinates, center);
            FitImgGreen = gaussFun2D_cov_fixedC(parFitGreen, coordinates);
            AA_OUT_COEF(5,i) = parFitGreen(2);
            AA_OUT_COEF(6,i) = (parFitGreen(1)*parFitGreen(2)*sqrt(2*pi))^2;
        end
        
    end
    
    
    
    % Plot Stuff
    scrsz = get(0,'ScreenSize');
    figure('Position',[1 scrsz(4)/2 2*scrsz(3)/3 scrsz(4)])
    
    % Plot first channel and track
    hRed = subplot(2,2,1);
    subimage(ImRed);
    title('First Channel');
    hold on;
    plot(tracks_corr{i}(:,2)/pixelSize,tracks_corr{i}(:,3)/pixelSize,'r');
    hold off;
    
    % Plot second channel and region for particle size estimation
    hGreen = subplot(2,2,2);
    if nChannels == 2
        subimage(ImGreen);
        title('Second Channel');
    else
        subimage(ImRed);
        title('First Channel');
    end
    
    if MeasureSize
        hold on;
        rectangle('Position',[xsub(1), ysub(1), fitWindow,fitWindow], 'EdgeColor', 'y');
        hold off;
    end
    
    
    % Plot MSD
    subplot(2,2,3)
    errorbar(AA_tlist, AA_OUT_MSD(:,i), AA_OUT_MSDerr(:,i),'ok');
    title({'MSD Plot';['D = ', num2str(AA_OUT_COEF(1,i)), 'µm^2/s'];...
        ['\alpha = ',num2str(AA_OUT_COEF(2,i))]});
    hold on;
    plot(AA_tlist, MSD(:,4),'r');
    hold off
    legend('Data', 'AD Fit');
    xlabel('Time [s]');
    ylabel('MSD [µm^2]');
    
    
    % Plot horizontal and vertical profiles for fitted particle
    if MeasureSize
        subplot(4,2,6)
        
        
        plot(ImRedsub(fitWindow/2+1,:),'or');
        hold on;
        plot(FitImgRed(fitWindow/2+1,:),'k');
        
        if nChannels == 2
            plot(ImGreensub(fitWindow/2+1,:),'og');
            
            plot(FitImgGreen(fitWindow/2+1,:),'k');
            
            title({'Horizontal profile of aggregate';...
                ['\sigma_r = ', num2str(AA_OUT_COEF(3,i)),...
                '; \sigma_g = ', num2str(AA_OUT_COEF(5,i))]});
            
        else
            title({'Horizontal profile of aggregate';...
                ['\sigma_r = ', num2str(AA_OUT_COEF(3,i))]});
        end
        hold off;
        
        subplot(4,2,8)
        
        plot(ImRedsub(:,fitWindow/2+1),'or');
        hold on;
        plot(FitImgRed(:,fitWindow/2+1),'k');
        if nChannels == 2
            
            plot(ImGreensub(:,fitWindow/2+1),'og');
            plot(FitImgGreen(:,fitWindow/2+1),'k');
            
            title({'Vertical profile of aggregate';...
                ['I_r = ', num2str(AA_OUT_COEF(4,i)),...
                '; I_g = ', num2str(AA_OUT_COEF(6,i))]});
        else
            
            
            title({'Vertical profile of aggregate';...
                ['I_r = ', num2str(AA_OUT_COEF(4,i))]});
        end
        
        hold off;
        
        
        
    end
    
    answer= questdlg('Next Track?');
    
    switch answer
        case 'Yes'
            close all;
        case 'No'
            break;
        case 'Cancel'
            break;
    end
end






