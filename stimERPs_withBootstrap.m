function stimERPs_withBootstrap(PatientID)
    
    globalemotionpath = '/Users/Aria/Documents/Emotion';
    patientdatapath = fullfile(globalemotionpath,'Data',PatientID);
    savedir = fullfile(patientdatapath,'AnalyzedData','SufficiencyStimMaps');

    computeMaps = 1;
    plotIndividualTimecourses = 0;
    fs = 100;


    if computeMaps
        %Load sufficiency stim struct
        [file, path] = uigetfile(patientdatapath);
        load(fullfile(path,file),'out_movie');
        
    %     %throw out sham stims
        out_movie = out_movie(logical([out_movie.stimApplied]));
        
    
       
    
    %     %temp for ec286 b16:
    %     goodinds = logical(ones(96,1));
    %     goodinds([19 33 43 44 47 50 51 53 57 64 69 79 89 94]) = logical(0);
    % %     out_movie = out_movie(goodinds);
        blocks = [out_movie.block];
        % out_movie = out_movie(blocks == 14);
        
        %Get indices that will contribute to each time series
        stimInfo = [out_movie.movie_info];
    % 
    %     target_category = extractfield(stimInfo,'pilot_emotion');
    %     out_movie = out_movie(find(strcmp(target_category,'Disgust')));
    
    %     stimInfo = [out_movie.movie_info];
    
        sources = {stimInfo(:).source};
        targets = {stimInfo(:).target};
        % sources = extractfield(stimInfo,'source');
        % targets = extractfield(stimInfo,'target');
        frequencies = [stimInfo.frequency];
        currents = [stimInfo.current];
    
        unique_sources = unique(sources);
    
        allMapInfo = struct();
        %Number of permutations for bootstrap
        numFolds = 1000;
        numMaps = 1;
        current_threshold = 1;
        metric = 'phasic';
        maxCurrent = 10;
        minCurrent = 0;
        maxDuration = 20;
        for i = 1:length(unique_sources)
            thisSource_inds = find(strcmp(sources,unique_sources{i}));
            thisSource_targets = unique(targets(thisSource_inds),'stable');
            for j = 1:length(thisSource_targets)
                thisSourceTargetPair_inds = intersect(find(strcmp(sources,unique_sources{i})),find(strcmp(targets,thisSource_targets{j})),'stable');
    
                %There may be 1Hz pulse trains across this same pair.
                thisSourceTargetPair_freqs = unique(frequencies(thisSourceTargetPair_inds));
    %             thisSourceTargetPair_freqs = thisSourceTargetPair_freqs(thisSourceTargetPair_freqs > 1);
                
                %Make an ERP for each source/target pair and frequency
                for k = 1:length(thisSourceTargetPair_freqs)
                    thisSourceTargetFreq_inds = intersect(thisSourceTargetPair_inds,find(frequencies == thisSourceTargetPair_freqs(k)),'stable');
                    allMapInfo(numMaps).inds = thisSourceTargetFreq_inds;
                    allMapInfo(numMaps).currents = currents(thisSourceTargetFreq_inds);
                    allMapInfo(numMaps).source = unique_sources{i};
                    allMapInfo(numMaps).target = thisSource_targets{j};
                    allMapInfo(numMaps).frequency = thisSourceTargetPair_freqs(k);
                    timecourseLengths = [];
                    for l = 1:length(thisSourceTargetFreq_inds)
                        if(current_threshold && ((out_movie(thisSourceTargetFreq_inds(l)).movie_info.current > maxCurrent) || (out_movie(thisSourceTargetFreq_inds(l)).movie_info.current < minCurrent)))
                            continue;
                        end
                        switch metric
                            case 'phasic'
                                timecourseLengths = [timecourseLengths length(out_movie(thisSourceTargetFreq_inds(l)).processed_biopacdata.phasic.rawData)];
                            case 'rawEDA'
                                timecourseLengths = [timecourseLengths length(out_movie(thisSourceTargetFreq_inds(l)).biopacdata(6,:))];
                            case 'tonic'
                                timecourseLengths = [timecourseLengths length(out_movie(thisSourceTargetFreq_inds(l)).processed_biopacdata.tonic.rawData)];
                            case 'RSA'
                                timecourseLengths = [timecourseLengths length(out_movie(thisSourceTargetFreq_inds(l)).processed_biopacdata.RSA.RSA_melSpec_12s.rawData)];
                            case 'IBI'
                                timecourseLengths = [timecourseLengths length(out_movie(thisSourceTargetFreq_inds(l)).processed_biopacdata.IBI.pchip.rawData)];
                        end
                    end
    
                    if(isempty(timecourseLengths))
                        continue;
                    end
    
                    minTimecourseLength = min(timecourseLengths);
                    minTimecourseLength = min(minTimecourseLength,maxDuration*fs);
    
                    allMapInfo(numMaps).timecourses = zeros(length(timecourseLengths),minTimecourseLength);
                    currTimecourse = 1;
                    for l = 1:length(thisSourceTargetFreq_inds)
                        if(current_threshold && ((out_movie(thisSourceTargetFreq_inds(l)).movie_info.current > maxCurrent) || (out_movie(thisSourceTargetFreq_inds(l)).movie_info.current < minCurrent)))
                            continue
                        end
                        % if(currTimecourse == 2)
                        %     currTimecourse = currTimecourse + 1;
                        %     continue;
                        % end
                        switch metric
                            case 'phasic'
                                allMapInfo(numMaps).timecourses(currTimecourse,:) = out_movie(thisSourceTargetFreq_inds(l)).processed_biopacdata.phasic.rawData(1:minTimecourseLength);
                            case 'rawEDA'
                                allMapInfo(numMaps).timecourses(currTimecourse,:) = out_movie(thisSourceTargetFreq_inds(l)).biopacdata(6,1:minTimecourseLength);
                            case 'tonic'
                                allMapInfo(numMaps).timecourses(currTimecourse,:) = out_movie(thisSourceTargetFreq_inds(l)).processed_biopacdata.tonic.rawData(1:minTimecourseLength);                            
                            case 'RSA'
                                allMapInfo(numMaps).timecourses(currTimecourse,:) = out_movie(thisSourceTargetFreq_inds(l)).processed_biopacdata.RSA.RSA_melSpec_12s.rawData(1:minTimecourseLength);  
                            case 'IBI'
                                allMapInfo(numMaps).timecourses(currTimecourse,:) = 60./out_movie(thisSourceTargetFreq_inds(l)).processed_biopacdata.IBI.pchip.rawData(1:minTimecourseLength);
                        end
                        allMapInfo(numMaps).timecourses(currTimecourse,:) = allMapInfo(numMaps).timecourses(currTimecourse,:) - mean(allMapInfo(numMaps).timecourses(currTimecourse,1:3*fs));
                        currTimecourse = currTimecourse + 1;
                    end
    
                    allMapInfo(numMaps).meanTimecourse = mean(allMapInfo(numMaps).timecourses,1);
                    allMapInfo(numMaps).semTimecourse = std(allMapInfo(numMaps).timecourses,1)/sqrt(size(allMapInfo(numMaps).timecourses,1));
                    [allMapInfo(numMaps).controlTimecourse, allMapInfo(numMaps).controlSemTimecourse] = stimBootstrap(allMapInfo(numMaps).timecourses,numFolds,fs);
    
                    numMaps = numMaps + 1;
                end
    
            end
        end
    
        mkdir(savedir);
        % save(fullfile(savedir,[PatientID '_B' num2str(out_movie(1).block) '_SufficiencyStim_PhasicMap']),'allMapInfo','-v7.3');
    end

    % load(fullfile(savedir,[PatientID '_B22_SufficiencyStim_PhasicMap']),'allMapInfo');
    % load(fullfile(savedir,[PatientID '_B' num2str(out_movie(1).block) '_SufficiencyStim_PhasicMap']),'allMapInfo');

    close all;
    bigfigN(1)
    tiledlayout('flow');
    addpath('/Users/Aria/Documents/Emotion/Code/');
    colors = {'green','red'};
    for i = 1:length(allMapInfo)
        nexttile;
% 
        p_real = stdshade_MeanSem(allMapInfo(i).meanTimecourse,allMapInfo(i).semTimecourse,.3,'red');
%         p_real = stdshade_MeanSem(allMapInfo(i).meanTimecourse,allMapInfo(i).semTimecourse,.3,colors{i});
        hold on;
        p_control = stdshade_MeanSem(allMapInfo(i).controlTimecourse,allMapInfo(i).controlSemTimecourse,.1,'black');
        if(plotIndividualTimecourses)
            for j = 1:size(allMapInfo(i).timecourses,1)
                plot(allMapInfo(i).timecourses(j,:),colors{i});
            end
        end
        title([allMapInfo(i).source ' to ' allMapInfo(i).target ' (' num2str(size(allMapInfo(i).timecourses,1)) ' stims)']);
        ylim([-.1 .4]);
        % ylim([-4 5]);
        xline(3*fs,'-','Stim Applied');
        xticks(0:fs:size(allMapInfo(i).timecourses,2));
        xticklabels(0:floor(size(allMapInfo(i).timecourses,2)/fs));
        xlabel('Time (s)');
        ylabel('Raw EDA change (μS)');

 
    end

%     saveas(gcf,fullfile(savedir,[PatientID '_B' num2str(out_movie(1).block) '_SufficiencyStim_' metric '_uncentered.png']));

end

function [meanTimecourse,stdTimecourse] = stimBootstrap(inputTimecourses,numFolds,fs)

    allPermuted = zeros(numFolds,length(inputTimecourses));

%     %If the number of folds exceeds the length of the timecourse, we'd get
%     %pigeonholed into duplicate permutations.
%     numFolds = min(numFolds,length(inputTimecourse));

    for i = 1:numFolds
        cycleDistances = randi(size(inputTimecourses,2),size(inputTimecourses,1));
        thisFoldShiftedTimecourses = zeros(size(inputTimecourses,1),size(inputTimecourses,2));
        for j = 1:size(thisFoldShiftedTimecourses,1)
            thisFoldShiftedtimecourses(j,:) = circshift(inputTimecourses(j,:),cycleDistances(j));
            thisFoldShiftedtimecourses(j,:) = thisFoldShiftedtimecourses(j,:) - mean(thisFoldShiftedtimecourses(j,1:3*fs));
        end
        
        allPermuted(i,:) = mean(thisFoldShiftedtimecourses,1);

    end

    [stdTimecourse meanTimecourse] = std(allPermuted,[],1);
    semTimecourse = stdTimecourse/sqrt(numFolds);
    
end

function bigfigN(N,h,w)
%opens a figure window of the specified size. If (height) h= 1 and (width) w=1 then the figure is the size of
%computer screen screen.
%N = figure number
%h = fraction of total screen height.
%w = fraction of tital screen width. 

%Hullett

if nargin == 1
    h = 1;
    w = 1;
end

scrsz = get(0,'ScreenSize');
fign = figure(N);
set(fign,'Position',[1 1 w*round(scrsz(3)) h*round(scrsz(4))]);
end