function stimERPs_withBootstrap_v3(PatientID)
    %Changes from v2:
    % - Does not include descending electrode # stims (e.g. OFC3-2), unless
    % the majority of stims were done that way. Uses either all ascending
    % pairs or all descending pairs.
    %TODO: Set metrics below: metric, yMin (for plotting), yMax (for
    %plotting)
    metric = 'phasic'
    yMin = -2;
    yMax = 6;
    Nboot = 3000;

    globalemotionpath = '/Users/Aria/Documents/Emotion';
    patientdatapath = fullfile(globalemotionpath,'Data',PatientID);
    savedir = fullfile(patientdatapath,'AnalyzedData','SufficiencyStimMaps','v3');

    fs = 100;
    stim_sites = {'AC','PC','OFCA','OFCP','A','AH','PH','AI','PI'};
    %Choose sufficiency stim struct
    [file, path] = uigetfile(patientdatapath);
    % load(fullfile(patientdatapath,'movie_out_structs','Sufficiency Stim','out_movie_100Hz_FINAL.mat'),'out_movie'); %For cleaned RSA data
    load(fullfile(path,file),'out_movie');
    
    %throw out sham stims
    out_movie = out_movie(logical([out_movie.stimApplied]));
    
    full_eda = [];
    full_rsa = [];

    blocks = [out_movie.block];
    % out_movie = out_movie(blocks == 14); %For EC286, outstruct includes additional non-sufficiency stim events

    %Gets full-length time series for null distribution calculation, including in between stim events
    switch metric
        case 'phasic'
            load(fullfile(patientdatapath,'Processed_BioPacData',strcat(PatientID,'_cvxEDA_SuffStim_FullTimeCourse.mat')),'cvxOutputFullTimeCourse');
            full_eda = cvxOutputFullTimeCourse.r';
        case 'RSA'
            load(fullfile(patientdatapath,'Processed_BioPacData',strcat(PatientID,'_Suff_RSA_Cleaned_FullTimeCourse_zScored.mat')),'rsa_cleaned_melSpec_8s_zData'); %For cleaned RSA data
            full_rsa = rsa_cleaned_melSpec_8s_zData;
            % full_rsa = rsa_cleaned_melSpec_16s_zData;
            % for i = 1:length(out_movie)
            %     l = min(length(out_movie(i).processed_biopacdata.RSA.RSA_melSpec_8s_cleaned.zData),2000);
            %     full_rsa = [full_rsa out_movie(i).processed_biopacdata.RSA.RSA_melSpec_8s_cleaned.zData(1:l)];
            % end
            % m = min(full_rsa)*-1
            % full_rsa = full_rsa + m;
            % for i = 1:length(out_movie)
            %     out_movie(i).processed_biopacdata.RSA.RSA_melSpec_8s_cleaned.zData = out_movie(i).processed_biopacdata.RSA.RSA_melSpec_8s_cleaned.zData + m;
            % end
    end

    %Get indices that will contribute to each time series
    stimInfo = [out_movie.movie_info];
    sources = {stimInfo(:).source};
    targets = {stimInfo(:).target};
    frequencies = [stimInfo.frequency];
    currents = [stimInfo.current];

    allMapInfo = struct();

    maxDuration = 20; %Maximum duration of stim event currently at 20s or length of shortest stim event otherwise
    for i = 1:7 %Sort by current
        numBins = 1;
        minCurrent = i - 1;
        maxCurrent = i;
        if i == 7
            thisCurrent_inds = find(currents > minCurrent); %Final map is 6mA and above
        else
            thisCurrent_inds = find(currents > minCurrent & currents <= maxCurrent); %Current range is (i-1,i]
        end
        allMapInfo(i).current = i; 
        allMapInfo(i).inds = thisCurrent_inds; 
        allMapInfo(i).bins = struct(); %Each bin is a unique source-target-frequency combination for this current range
        thisCurrent_sources = unique(sources(thisCurrent_inds),'stable'); %Doesn't matter if this is in any particular order. Will create a bins trait that stores source name
        if isempty(thisCurrent_sources)
            allMapInfo(i).hasData = 0;
        else
            allMapInfo(i).hasData = 1;
        end
        for j = 1:length(thisCurrent_sources) %Sort by source electrode
            source_elec = char(thisCurrent_sources(j));
            thisCurrentSource_inds = intersect(thisCurrent_inds,find(strcmp(sources,thisCurrent_sources(j)))); %default sorted order
            thisCurrentSource_targets = unique(targets(thisCurrentSource_inds),'stable'); %should mostly give the adjacent +1 electrode, flag cases that are not so
            for k = 1:length(thisCurrentSource_targets) %Sort by target electrode
                target_elec = char(thisCurrentSource_targets(k));
                thisCurrentSourceTarget_inds = intersect(thisCurrentSource_inds,find(strcmp(targets,thisCurrentSource_targets(k))),'stable'); 
                thisCurrentSourceTarget_frequencies = unique(frequencies(thisCurrentSourceTarget_inds),'stable');
                thisCurrentSourceTarget_frequencies = thisCurrentSourceTarget_frequencies(thisCurrentSourceTarget_frequencies > 1); %Remove 1Hz events
                if isempty(thisCurrentSourceTarget_frequencies) %Skips this bin for plotting if there were ONLY 1Hz events at this source-target pair
                    allMapInfo(i).hasData = 0;
                end
                for l = 1:length(thisCurrentSourceTarget_frequencies) %Sort by frequency
                    thisCurrentSourceTargetFreq_inds = intersect(thisCurrentSourceTarget_inds,find(frequencies==thisCurrentSourceTarget_frequencies(l)));
                    [flag source_elec target_elec] = checkAdjacency(source_elec,target_elec); %Checks if source and target electrodes are directly adjacent in ascending order (default)
                    % [flag source_elec target_elec] = checkAdjacency(source_elec,target_elec,0); %Use if default is descending electrode order
                    if flag || thisCurrentSourceTarget_frequencies(l) ~= 100
                        allMapInfo(i).bins(numBins).flagged = 1; %Flag this stim event if non-adjacent electrodes or if frequency setting is different
                    else
                        allMapInfo(i).bins(numBins).flagged = 0;
                    end
                    s_ind = getElecInd(source_elec);
                    t_ind = getElecInd(target_elec);
                    allMapInfo(i).bins(numBins).frequency = thisCurrentSourceTarget_frequencies(l);
                    allMapInfo(i).bins(numBins).source = source_elec;
                    allMapInfo(i).bins(numBins).target = target_elec; 
                    allMapInfo(i).bins(numBins).numStims = length(thisCurrentSourceTargetFreq_inds);
                    allMapInfo(i).bins(numBins).sourceLocation = checkAnatomical(source_elec(2:s_ind-1),stim_sites); %Corrects electrode site to fit standard names
                    allMapInfo(i).bins(numBins).targetLocation = checkAnatomical(target_elec(2:t_ind-1),stim_sites);
                    allMapInfo(i).bins(numBins).sourcePosition = str2num(source_elec(s_ind:end));
                    allMapInfo(i).bins(numBins).targetPosition = str2num(target_elec(t_ind:end));
                    if ~ismember(allMapInfo(i).bins(numBins).sourceLocation,stim_sites)
                        allMapInfo(i).bins(numBins).flagged = 1; %Flag abnormal sites
                    end
                    timecourseLengths = [];
                    for m = 1:length(thisCurrentSourceTargetFreq_inds)
                        switch metric
                            case 'phasic'
                                timecourseLengths = [timecourseLengths length(out_movie(thisCurrentSourceTargetFreq_inds(m)).processed_biopacdata.phasic.rawData)];
                            case 'rawEDA'
                                timecourseLengths = [timecourseLengths length(out_movie(thisCurrentSourceTargetFreq_inds(m)).biopacdata(6,:))];
                            case 'tonic'
                                timecourseLengths = [timecourseLengths length(out_movie(thisCurrentSourceTargetFreq_inds(m)).processed_biopacdata.tonic.rawData)];
                            case 'RSA'
                                timecourseLengths = [timecourseLengths length(out_movie(thisCurrentSourceTargetFreq_inds(m)).processed_biopacdata.RSA.RSA_melSpec_8s.rawData)];
                                % timecourseLengths = [timecourseLengths length(out_movie(thisCurrentSourceTargetFreq_inds(m)).processed_biopacdata.RSA.RSA_melSpec_8s_cleaned.zData)];
                            case 'IBI'
                                timecourseLengths = [timecourseLengths length(out_movie(thisCurrentSourceTargetFreq_inds(m)).processed_biopacdata.IBI.pchip.rawData)];
                        end
                    end
                    if(isempty(timecourseLengths)) %this happens if there are no indices - is that possible? Do I need this?
                        continue;
                    end
                    minTimecourseLength = min(timecourseLengths); 
                    minTimecourseLength = min(minTimecourseLength,maxDuration*fs);
                    timecourses = zeros(length(thisCurrentSourceTargetFreq_inds),minTimecourseLength);
                    for m = 1:length(thisCurrentSourceTargetFreq_inds)
                        switch metric
                            case 'phasic'
                                timecourses(m,:) = out_movie(thisCurrentSourceTargetFreq_inds(m)).processed_biopacdata.phasic.rawData(1:minTimecourseLength);
                            case 'rawEDA'
                                timecourses(m,:) = out_movie(thisCurrentSourceTargetFreq_inds(m)).biopacdata(6,1:minTimecourseLength);
                            case 'tonic'
                                timecourses(m,:) = out_movie(thisCurrentSourceTargetFreq_inds(m)).processed_biopacdata.tonic.rawData(1:minTimecourseLength);                            
                            case 'RSA'
                                timecourses(m,:) = out_movie(thisCurrentSourceTargetFreq_inds(m)).processed_biopacdata.RSA.RSA_melSpec_8s_cleaned.zData(1:minTimecourseLength)'; %To use cleaned data
                                % timecourses(m,:) = out_movie(thisCurrentSourceTargetFreq_inds(m)).processed_biopacdata.RSA.RSA_melSpec_16s_cleaned.zData(1:minTimecourseLength)'; %To use cleaned data
                            case 'IBI' 
                                timecourses(m,:) = 60./out_movie(thisCurrentSourceTargetFreq_inds(m)).processed_biopacdata.IBI.pchip.rawData(1:minTimecourseLength);
                        end  
                        timecourses(m,:) = timecourses(m,:) - mean(timecourses(m,1:3*fs));
                    end
                    allMapInfo(i).bins(numBins).rawTimecourses = timecourses;
                    allMapInfo(i).bins(numBins).meanTimecourse = mean(timecourses,1);
                    allMapInfo(i).bins(numBins).semTimecourse = std(timecourses,1)/sqrt(length(thisCurrentSourceTargetFreq_inds));
                    switch metric
                        case 'phasic'
                            % [allMapInfo(i).bins(numBins).controlTimecourse allMapInfo(i).bins(numBins).controlSemTimecourse] = stimBootstrap(full_eda,length(thisCurrentSourceTargetFreq_inds),fs,minTimecourseLength);
                            [allMapInfo(i).bins(numBins).controlTimecourse allMapInfo(i).bins(numBins).controlSemTimecourse] = stimBootstrap_v2(full_eda,length(thisCurrentSourceTargetFreq_inds),fs,minTimecourseLength,Nboot);
                        case 'RSA'
                            % [allMapInfo(i).bins(numBins).controlTimecourse allMapInfo(i).bins(numBins).controlSemTimecourse] = stimBootstrap(full_rsa,length(thisCurrentSourceTargetFreq_inds),fs,minTimecourseLength);
                            [allMapInfo(i).bins(numBins).controlTimecourse allMapInfo(i).bins(numBins).controlSemTimecourse] = stimBootstrap_v2(full_rsa,length(thisCurrentSourceTargetFreq_inds),fs,minTimecourseLength,Nboot);
                    end  

                    %Cross-checking output against true bootstrapping:
                    out_movie_thisBin = out_movie(thisCurrentSourceTargetFreq_inds);
                    switch metric
                        case 'phasic'
                            [sitemn mn_boot std_boot sitemnboot t std_boot_pos std_boot_neg] = get_stim_mean_nulldist_EDA(out_movie_thisBin, out_movie,Nboot);
                        case 'RSA'
                            [sitemn mn_boot std_boot sitemnboot t std_boot_pos std_boot_neg] = get_stim_mean_nulldist_RSA(out_movie_thisBin, out_movie,Nboot);
                    end
                    allMapInfo(i).bins(numBins).sitemn = sitemn;
                    allMapInfo(i).bins(numBins).mn_boot = mn_boot;
                    allMapInfo(i).bins(numBins).std_boot_pos = std_boot_pos;
                    allMapInfo(i).bins(numBins).std_boot_neg = std_boot_neg;
                    allMapInfo(i).bins(numBins).sitemnboot = sitemnboot;
                    allMapInfo(i).bins(numBins).t = t;

                    numBins = numBins + 1;
                end
            end
        end
    end
    
    mkdir(savedir);
    close all;
    rows = dictionary(stim_sites,[1:9]);
    elec_paths = {'9','8','7','6','5','4','3','2','1'}; %9-10, 8-9, ... , 1-2 (L), (R) 1-2, 2-3, ... , 8-9, 9-10
    columns = dictionary(str2double(elec_paths),[1:9]);

    addpath('/Users/Aria/Documents/Emotion/Code/');
    colors = jet(length(allMapInfo));
    fn = 1;
    for i = 1:length(allMapInfo)
        if ~allMapInfo(i).hasData %Skip this map if no stim events occurred in this current range
            continue;
        end

        %Figure set up - same for every page
        figi = bigfigN(fn)
        tiledlayout(9,18,"TileSpacing",'tight',"Padding",'tight',"OuterPosition",[.03 0 .97 .97]);
        for j = 1:length(stim_sites)
            dim = [.005 (9-j)*0.1 0.1 0.1];
            str = stim_sites(j);
            annotation('textbox',dim,'String',str,'EdgeColor','none','FontSize',8)
        end
        for j = 1:length(elec_paths)
            str = strcat(elec_paths(j), '-', num2str(str2double(elec_paths(j)) + 1));
            dim = [0.04+(0.05441*(j-1)) 0.9 0.1 0.1];
            annotation('textbox',dim,'String',str,'EdgeColor','none','FontSize',9)
            dim = [0.965-(0.05441*(j-1)) 0.9 0.1 0.1];
            annotation('textbox',dim,'String',str,'EdgeColor','none','FontSize',9)
        end

        flagged = [];
        for j = 1:length(allMapInfo(i).bins)
            if allMapInfo(i).bins(j).flagged %Skip flagged bins, plot on separate page immediately after this page
                flagged = [flagged j];
                continue;
            end

            %Find correct row and column for this bin:
            laterality = allMapInfo(i).bins(j).source(1);
            if laterality == 'L'
                laterality = 0;
            elseif laterality == 'R'
                laterality = 1;
            end
            row = rows(cellstr(allMapInfo(i).bins(j).sourceLocation));
            col = columns(allMapInfo(i).bins(j).sourcePosition) + laterality * (2 * allMapInfo(i).bins(j).sourcePosition - 1);       
            tilenum = (row - 1) * 18 + col;
            nexttile(tilenum);


            %Plot timecourses mean and sem:
            p_real = stdshade_MeanSem(allMapInfo(i).bins(j).meanTimecourse,0,.3,colors(i,:));
            hold on;
            p_control = stdshade_MeanSem(allMapInfo(i).bins(j).controlTimecourse,allMapInfo(i).bins(j).controlSemTimecourse,.1,'black');

            %OR plot individual timecourses, no sem:
            % hold on;
            % for k = 1:size(allMapInfo(i).bins(j).rawTimecourses,1)
            %     plot(allMapInfo(i).bins(j).rawTimecourses(k,:),'Color',colors(i,:));
            %     hold on;
            % end

            title([num2str(allMapInfo(i).bins(j).numStims) ' stim(s)']);
            xline(3*fs,'-','LineWidth',2) %Stim Applied timestamp
            ylim([yMin yMax]);
            a = get(gca,'XTickLabel');
            set(gca,'XTickLabel',a,'Fontsize',5);
            a = get(gca,'YTickLabel');
            set(gca,'YTickLabel',a,'Fontsize',5);
        end
        figname = [string(i-1) 'to' string(i) 'mA'];
        saveas(figi, strrep(strjoin([{patientdatapath} '/AnalyzedData/SufficiencyStimMaps/v3/' PatientID '_SuffStim_' metric '_' figname]), ' ', ''), 'png');
        % saveas(figi, strrep(strjoin([{patientdatapath} '/AnalyzedData/SufficiencyStimMaps/' PatientID '_SuffStim_' metric '_' figname '_notAveraged']), ' ', ''), 'png');
        fn = fn + 1;
        
        %Return to flagged bins and plot on new page, no specific format
        if ~isempty(flagged)
            flagi = bigfigN(fn)
            tiledlayout('flow');
            for k = 1:length(flagged)
                nexttile;
                p_real = stdshade_MeanSem(allMapInfo(i).bins(flagged(k)).meanTimecourse,allMapInfo(i).bins(flagged(k)).semTimecourse,.3,colors(i,:));
                hold on;
                p_control = stdshade_MeanSem(allMapInfo(i).bins(flagged(k)).controlTimecourse,allMapInfo(i).bins(flagged(k)).controlSemTimecourse,.1,'black');
                title([allMapInfo(i).bins(flagged(k)).source ' to ' allMapInfo(i).bins(flagged(k)).target ' (' num2str(allMapInfo(i).bins(flagged(k)).numStims) ' stim(s), ' num2str(allMapInfo(i).bins(flagged(k)).frequency) ' Hz)']);
                xline(3*fs,'-','Stim Applied');
                ylim([yMin yMax]);
                xticks(0:fs:size(allMapInfo(i).bins(flagged(k)).meanTimecourse,2));
                xticklabels(0:floor(size(allMapInfo(i).bins(flagged(k)).meanTimecourse,2)/fs));
                xlabel('Time (s)');
            end
            fn = fn + 1;
            figname = [string(i-1) 'to' string(i) 'mA_Wildcards'];
            saveas(flagi, strrep(strjoin([{patientdatapath} '/AnalyzedData/SufficiencyStimMaps/v3/' PatientID '_SuffStim_' metric '_' figname]), ' ', ''), 'png');
            % saveas(flagi, strrep(strjoin([{patientdatapath} '/AnalyzedData/SufficiencyStimMaps/' PatientID '_SuffStim_' metric '_' figname '_smallY']), ' ', ''), 'png');
        end

        %Cross-checking output against true bootstrapping:
        figi = bigfigN(fn)
        tiledlayout(9,18,"TileSpacing",'tight',"Padding",'tight',"OuterPosition",[.03 0 .97 .97]);
        for j = 1:length(stim_sites)
            dim = [.005 (9-j)*0.1 0.1 0.1];
            str = stim_sites(j);
            annotation('textbox',dim,'String',str,'EdgeColor','none','FontSize',8)
        end
        for j = 1:length(elec_paths)
            str = strcat(elec_paths(j), '-', num2str(str2double(elec_paths(j)) + 1));
            dim = [0.04+(0.05441*(j-1)) 0.9 0.1 0.1];
            annotation('textbox',dim,'String',str,'EdgeColor','none','FontSize',9)
            dim = [0.965-(0.05441*(j-1)) 0.9 0.1 0.1];
            annotation('textbox',dim,'String',str,'EdgeColor','none','FontSize',9)
        end
        for j = 1:length(allMapInfo(i).bins)
            if allMapInfo(i).bins(j).flagged || isempty(allMapInfo(i).bins(j).sitemn)%Skip flagged bins
                % %For EC286 only:
                % figi = bigfigN(fn) %Reset plot
                % plot(allMapInfo(i).bins(j).t,mean(allMapInfo(i).bins(j).sitemn,1),'LineWidth',1.5,'Color',colors(i,:));
                % hold on;
                % stdshade_mn_sem_computed_shadeOnly(allMapInfo(i).bins(j).t,allMapInfo(i).bins(j).mn_boot,allMapInfo(i).bins(j).std_boot_pos,allMapInfo(i).bins(j).std_boot_neg,[.2 .2 .2],.1,.001)
                % title([allMapInfo(i).bins(flagged(k)).source ' to ' allMapInfo(i).bins(flagged(k)).target ' (' num2str(allMapInfo(i).bins(flagged(k)).numStims) ' stim(s), ' num2str(allMapInfo(i).bins(flagged(k)).frequency) ' Hz)']);
                % xline(3,'-','Stim Applied');
                % ylim([yMin yMax]);
                % xlabel('Time (s)'); 
                % fn = fn + 1;
                % figname = [string(i-1) 'to' string(i) 'mA_Wildcards'];
                % saveas(flagi, strrep(strjoin([{patientdatapath} '/AnalyzedData/SufficiencyStimMaps/v3/' PatientID '_SuffStim_' metric '_' figname '_PH']), ' ', ''), 'png');

                continue;
            end

            %Find correct row and column for this bin:
            laterality = allMapInfo(i).bins(j).source(1);
            if laterality == 'L'
                laterality = 0;
            elseif laterality == 'R'
                laterality = 1;
            end
            row = rows(cellstr(allMapInfo(i).bins(j).sourceLocation));
            col = columns(allMapInfo(i).bins(j).sourcePosition) + laterality * (2 * allMapInfo(i).bins(j).sourcePosition - 1);       
            tilenum = (row - 1) * 18 + col;
            nexttile(tilenum);

            %Plot timecourses mean and sem:
            plot(allMapInfo(i).bins(j).t,mean(allMapInfo(i).bins(j).sitemn,1),'LineWidth',1.5,'Color',colors(i,:));
            hold on;
            stdshade_mn_sem_computed_shadeOnly(allMapInfo(i).bins(j).t,allMapInfo(i).bins(j).mn_boot,allMapInfo(i).bins(j).std_boot_pos,allMapInfo(i).bins(j).std_boot_neg,[.2 .2 .2],.1,.001)

            title([num2str(allMapInfo(i).bins(j).numStims) ' stim(s)']);
            xline(3,'-','LineWidth',2) %Stim Applied timestamp
            ylim([yMin yMax]);
            a = get(gca,'XTickLabel');
            set(gca,'XTickLabel',a,'Fontsize',5);
            a = get(gca,'YTickLabel');
            set(gca,'YTickLabel',a,'Fontsize',5);
        end
        figname = [string(i-1) 'to' string(i) 'mA'];
        saveas(figi, strrep(strjoin([{patientdatapath} '/AnalyzedData/SufficiencyStimMaps/v3/' PatientID '_SuffStim_' metric '_' figname '_PH']), ' ', ''), 'png');
        fn = fn + 1;
    end
end

function [meanTimecourse,stdTimecourse] = stimBootstrap(inputTimecourse,numFolds,fs,minTimeCourseLength)

    allPermuted = zeros(numFolds,minTimeCourseLength);

%     %If the number of folds exceeds the length of the timecourse, we'd get
%     %pigeonholed into duplicate permutations.
%     numFolds = min(numFolds,length(inputTimecourse));
    for i = 1:numFolds
        thisFoldShiftedtimecourse = circshift(inputTimecourse,randi(length(inputTimecourse)));
        thisFoldShiftedtimecourse = thisFoldShiftedtimecourse(1:minTimeCourseLength);
        thisFoldShiftedtimecourse = thisFoldShiftedtimecourse - mean(thisFoldShiftedtimecourse(1:3*fs));
        allPermuted(i,:) = thisFoldShiftedtimecourse;
    end

    [stdTimecourse meanTimecourse] = std(allPermuted,[],1); %Error bar: std * 2 +/- mean of means
    semTimecourse = stdTimecourse/sqrt(numFolds);
    
end

function [meanTimecourse,stdTimecourse] = stimBootstrap_v2(inputTimecourse,numEvents,fs,minTimeCourseLength,numFolds)
    
    allPermuted = zeros(numFolds,minTimeCourseLength);
%     %If the number of folds exceeds the length of the timecourse, we'd get
%     %pigeonholed into duplicate permutations.
%     numFolds = min(numFolds,length(inputTimecourse));
    for n = 1:numFolds
        sample_events = zeros(numEvents,minTimeCourseLength);
        for i = 1:numEvents
            rand_start = randi(length(inputTimecourse));
            if (rand_start + minTimeCourseLength) > length(inputTimecourse)
                l = length(inputTimecourse) - rand_start;
                thisFoldtimecourse = inputTimecourse(rand_start+1:end);
                thisFoldtimecourse = [thisFoldtimecourse inputTimecourse(1:minTimeCourseLength - l)];
            else
                thisFoldtimecourse = inputTimecourse(rand_start+1:rand_start+minTimeCourseLength);
            end
            % minTimeCourseLength
            % length(thisFoldtimecourse)
            thisFoldtimecourse = thisFoldtimecourse - mean(thisFoldtimecourse(1:3*fs));
            % thisFoldShiftedtimecourse = circshift(inputTimecourse,randi(length(inputTimecourse)));
            % thisFoldShiftedtimecourse = thisFoldShiftedtimecourse(1:minTimeCourseLength);
            % thisFoldShiftedtimecourse = thisFoldShiftedtimecourse - mean(thisFoldShiftedtimecourse(1:3*fs));
            sample_events(i,:) = thisFoldtimecourse;
        end
        [stdTimecourse meanTimecourse] = std(sample_events,[],1); 
        allPermuted(n,:) = meanTimecourse;
    end
    [stdTimecourse meanTimecourse] = std(allPermuted,[],1); %Error bar: meantimecourse2 +/- 2*stdtimecourse2
end

function fign = bigfigN(N,h,w)
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

function num = getElecInd(site)
    %Gets 1-10 depth electrode index based on electrode name
    if ~isempty(str2num(site(end-1)))
        num = length(site) - 1; %Returns index of electrode number
        % num = str2num(site(end-1:end)); %Returns electrode number
    else
        num = length(site);
        % num = str2num(site(end));
    end
end

function [flag,source_elec,target_elec] = checkAdjacency(source,target,ascending)
    %Checks if source and target are adjacent electrodes: same depth
    %electrode & adjacent 1-10 positions
    if nargin < 3
        ascending = 1;
    end
    source_elec = source;
    target_elec = target;
    source_num = getElecInd(source);
    target_num = getElecInd(target);
    if ~strcmp(source(1:source_num-1),target(1:target_num-1)) %Check same site first
        flag = 1;
        return;
    end
    if ascending == 1
        if str2num(source(source_num:end)) ~= str2num(target(target_num:end)) - 1 %E.g. 10->9
        % flag = 0;
            flag = 1;
        else
            flag = 0;
        end
        % source_elec = target; %As of 10/27/23, no longer treating 2-3 as equal to 3-2
        % target_elec = source;
    elseif str2num(source(source_num:end)) ~= str2num(target(target_num:end)) + 1
        flag = 1;
    else
        flag = 0;
        source_elec = target;
        target_elec = source;
    end
end
    
function site = checkAnatomical(site,stim_sites)
    %If no anterior/posterior orientation specified, assume anterior except
    %for OFC
    if ismember(site,stim_sites)
        return;
    elseif strcmp(site,'C') || strcmp(site,'CD') || strcmp(site,'ACD')
        site = 'AC';
    elseif strcmp(site,'PCD') || strcmp(site,'CM')
        site = 'PC';
    elseif strcmp(site,'OFC') || strcmp(site,'OFCD')
        site = 'OFCP';
    elseif strcmp(site,'H') || strcmp(site,'HD')
        site = 'AH';
    elseif strcmp(site,'I') || strcmp(site,'ID') || strcmp(site,'IA')
        site = 'AI';
    elseif strcmp(site,'IP')
        site = 'PI';
    end
end

    
