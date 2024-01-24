function plot_from_allMapInfo(PatientID)
    globalemotionpath = '/Users/Aria/Documents/Emotion';
    patientdatapath = fullfile(globalemotionpath,'Data',PatientID);

    %Manually select allMapInfo .mat file
    [file, path] = uigetfile(patientdatapath);
    load(fullfile(path,file),'allMapInfo');

    %TODO:
    yMin = -1;
    yMax = 4.5;
    m = [3 5 6]; %Choose which maps and corresponding bins you want to plot
    b = [2 5 1];
    plot_title = 'Mid-Cingulate Skin Conductance Response to Stimulation';
    y_title = 'Skin Conductance (z-score)';
    legend_labels = {'6 mA','5 mA','3 mA'};

    % colors = jet(length(allMapInfo));
    colors = hot(length(allMapInfo));
    colors(3,3) = 0.5000;
    fig1 = bigfigN(1)

    for i = 1:length(m)
        plot(allMapInfo(m(i)).bins(b(i)).t,allMapInfo(m(i)).bins(b(i)).sitemn,'Color',colors(i,:),'LineWidth',3);
        stdshade_mn_sem_computed_shadeOnly(allMapInfo(m(i)).bins(b(i)).t,allMapInfo(m(i)).bins(b(i)).mn_boot,allMapInfo(m(i)).bins(b(i)).std_boot_pos,allMapInfo(m(i)).bins(b(i)).std_boot_neg,colors(i,:),.03+(i*.01),.001);
        hold on;
        % L{i,1} = plot(nan,nan,'Color',colors(i,:));
        hold on;
    end
    xline(3,'-','Stim Applied');
    ylim([yMin yMax]);
    title(plot_title,'FontSize',16);
    xlabel('Time (s)','FontSize',18);
    ylabel(y_title,'FontSize',18);

    L1 = plot(nan,nan,'Color',colors(3,:),'LineWidth',3);
    L2 = plot(nan,nan,'Color',colors(2,:),'LineWidth',3);
    L3 = plot(nan,nan,'Color',colors(1,:),'LineWidth',3);

    legend([L1,L2,L3],legend_labels);
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