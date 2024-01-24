function stimERPs_withBootstrap_v4_mod(PatientID,metric,yMin,yMax,plotWindow,specWindow)
    %TRUE BOOTSTRAPPING ONLY

    %TODO: Set metrics below: metric, yMin (for plotting), yMax (for
    %plotting)
    % metric = 'phasic'
    % yMin = -3;
    % yMax = 3;
    Nboot = 3000;
    if nargin < 6
        specWindow = 0;
    end

    globalemotionpath = getenv('EMOTION_ROOT');
    patientdatapath = fullfile(globalemotionpath,'Data',PatientID);
    savedir = fullfile(patientdatapath,'AnalyzedData','SufficiencyStimMaps','v4');

    stim_sites = {'AC','PC','OFCA','OFCP','A','AH','PH','AI','PI'};
    %Choose sufficiency stim struct
    % [file, path] = uigetfile(patientdatapath);
    load(fullfile(patientdatapath,'movie_out_structs','Sufficiency Stim','out_movie_100Hz_withRSA_PPG.mat'),'out_movie');
    % load(fullfile(patientdatapath,'movie_out_structs','Sufficiency Stim','out_movie_100Hz_FINAL.mat'),'out_movie'); %For cleaned RSA data
    % load(fullfile(path,file),'out_movie');
    
    %throw out sham stims
    out_movie = out_movie(logical([out_movie.stimApplied]));
    
    blocks = [out_movie.block];

    %Get indices that will contribute to each time series
    stimInfo = [out_movie.movie_info];
    sources = {stimInfo(:).source};
    targets = {stimInfo(:).target};
    frequencies = [stimInfo.frequency];
    currents = [stimInfo.current];

    allMapInfo = struct();

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

                    out_movie_thisBin = out_movie(thisCurrentSourceTargetFreq_inds);
                    switch metric
                        case 'phasic'
                            [sitemn mn_boot std_boot sitemnboot t std_boot_pos std_boot_neg] = get_stim_mean_nulldist_EDA(out_movie_thisBin, out_movie,Nboot,plotWindow);
                        case 'RSA'
                            [sitemn mn_boot std_boot sitemnboot t std_boot_pos std_boot_neg] = get_stim_mean_nulldist_RSA_mod(out_movie_thisBin, out_movie,Nboot,plotWindow,specWindow);
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
    switch metric
        case 'phasic'
            save(fullfile(patientdatapath,'AnalyzedData',['stimERP_v4_' metric '_allMapInfo_' num2str(plotWindow) 'sec.mat']),'allMapInfo','-v7.3');
        case 'RSA'
            save(fullfile(patientdatapath,'AnalyzedData',['stimERP_v4_' metric '_allMapInfo_' num2str(plotWindow) '+' num2str(specWindow) 'sec.mat']),'allMapInfo','-v7.3');
    end
    
    mkdir(savedir);
    close all;
    rows = dictionary(stim_sites,[1:9]);
    elec_paths = {'9','8','7','6','5','4','3','2','1'}; %9-10, 8-9, ... , 1-2 (L), (R) 1-2, 2-3, ... , 8-9, 9-10
    columns = dictionary(str2double(elec_paths),[1:9]);

    addpath(genpath(fullfile(globalemotionpath,'Code')));
    colors = jet(length(allMapInfo));
    fn = 1;
    for i = 1:length(allMapInfo)
        if ~allMapInfo(i).hasData %Skip this map if no stim events occurred in this current range
            continue;
        end

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
                if strcmp(PatientID,'EC286')
                    figi = bigfigN(fn) %Reset plot
                    plot(allMapInfo(i).bins(j).t,mean(allMapInfo(i).bins(j).sitemn,1),'LineWidth',1.5,'Color',colors(i,:));
                    hold on;
                    stdshade_mn_sem_computed_shadeOnly(allMapInfo(i).bins(j).t,allMapInfo(i).bins(j).mn_boot,allMapInfo(i).bins(j).std_boot_pos,allMapInfo(i).bins(j).std_boot_neg,[.2 .2 .2],.1,.001)
                    title([allMapInfo(i).bins(j).source ' to ' allMapInfo(i).bins(j).target ' (' num2str(allMapInfo(i).bins(j).numStims) ' stim(s), ' num2str(allMapInfo(i).bins(j).frequency) ' Hz)']);
                    xline(3,'-','Stim Applied');
                    ylim([yMin yMax]);
                    xlabel('Time (s)'); 
                    fn = fn + 1;
                    figname = [string(i-1) 'to' string(i) 'mA_Wildcards'];
                    switch metric
                        case 'phasic'
                            saveas(figi, strrep(strjoin([{patientdatapath} '/AnalyzedData/SufficiencyStimMaps/v4/' PatientID '_SuffStim_' metric '_' figname '_' num2str(plotWindow) 'sec']), ' ', ''), 'png');
                        case 'RSA'
                            saveas(figi, strrep(strjoin([{patientdatapath} '/AnalyzedData/SufficiencyStimMaps/v4/' PatientID '_SuffStim_' metric '_' figname '_' num2str(plotWindow) '+' num2str(specWindow) 'sec']), ' ', ''), 'png');
                    end
                end
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
        switch metric
            case 'phasic'
                saveas(figi, strrep(strjoin([{patientdatapath} '/AnalyzedData/SufficiencyStimMaps/v4/' PatientID '_SuffStim_' metric '_' figname '_' num2str(plotWindow) 'sec']), ' ', ''), 'png');
            case 'RSA'
                saveas(figi, strrep(strjoin([{patientdatapath} '/AnalyzedData/SufficiencyStimMaps/v4/' PatientID '_SuffStim_' metric '_' figname '_' num2str(plotWindow) '+' num2str(specWindow) 'sec']), ' ', ''), 'png');
        end
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

    
