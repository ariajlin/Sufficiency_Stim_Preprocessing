function [sitemn mn_boot std_boot sitemnboot t std_boot_pos std_boot_neg] = get_stim_mean_nulldist_RSA(out_subSet, out_all,Nboot,indStop_sec,spec_window);

%% inputs
%out_all = outstruct for all the stimlation trials in a block or for a
    %participant. The null distribution is taken from this. 
%out_subSet = outstruct for just the stim trials at a specific site and
    %parameter range (for example: all the trials at OFC1-OFC1, 1-2mA, 100
    %Hz,5sec

%% outputs
% sitemn = the stimulus triggered mean for all the trials in out_subSet
% mn_boot = the mean of the null distribution (i.e. the mean of the set of
    % null means). 
% std_boot = the standard deviation of the null distribution (i.e. the standard deviation of the set of
    % null standard deviations);
% sitemnboot= a cell containing all the null destribution means (nthe
    % number is set by Nboot below). 
% t= time vector for tials 
% std_boot_pos = the true upper 95% confidence interval of the null distribution (i.e. approximated directly from  set of
    % null standard deviations);
% std_boot_neg = the true lower 95% confidence interval of the null distribution (i.e. approximated directly from  set of
    % null standard deviations)

FS = 100;%sample frequency
baseline = FS*3;%baseline to subtract. 
% indStop_sec = 17; %how many seconds to take off the end of movie.
% Nboot = 3000; %nuber of samples for the null distribution. 

%% first get the all the EDA (first ~20 sec of each stim trial) to draw from
%for the null distibution.k 
out = out_all;
 for j = 1:length(out)
%             if out(j).befaft(1) ~= out(j).processed_biopacdata.phasic.befaft(1)
%                 error('pre-movie baseline not the same')
%             end

        %resp =  out(j).processed_biopacdata.phasic.data;
        %resp = out(j).processed_biopacdata.phasic.rawData';
        switch spec_window
            case '8'
                resp = out(j).processed_biopacdata.RSA.RSA_melSpec_8s_cleaned.zData';
            case '16'
                resp = out(j).processed_biopacdata.RSA.RSA_melSpec_16s_cleaned.zData';
        end
        %resp = out(j).processed_biopacdata.phasic.zData'

        % ind_noise = find(resp> noise_rsa{pN,1});
        % resp(ind_noise) = 0;

        [rr cc] = size(resp);
        if rr > cc
            resp = resp';
        end

        
        out(j).resp = resp;  
 end
 
 
%get only ~20 seconds after stim
indrmv = [];
for j = 1:length(out)
    
    if length(out(j).resp) > indStop_sec*FS
        respj = out(j).resp(1:indStop_sec*FS);
        %respj = respj - mean(respj(1:baseline));
        out(j).resp = respj; 
    else
        indrmv = [indrmv j];
    end
 
%             
end

out(indrmv) = [];

 
%make minimum value 0 
resp_temp = [];
for j = 1:length(out)
resp_temp = [resp_temp out(j).resp];
end

%% now prosces the out_subSet in the same way. (set eda field = out.resp and truncate it to be ~20 seconds long). 
out = out_subSet;
for j = 1:length(out)

        switch spec_window
            case '8'
                resp = out(j).processed_biopacdata.RSA.RSA_melSpec_8s_cleaned.zData';
            case '16'
                resp = out(j).processed_biopacdata.RSA.RSA_melSpec_16s_cleaned.zData';
        end
        %resp = out(j).processed_biopacdata.phasic.zData'

        
        [rr cc] = size(resp);
        if rr > cc
            resp = resp';
        end

        
        out(j).resp = resp;  
 end
 
 
%get only ~20 seconds after stim
indrmv = [];
for j = 1:length(out)
    
    if length(out(j).resp) > indStop_sec*FS
        respj = out(j).resp(1:indStop_sec*FS);
        respj = respj - mean(respj(1:baseline));
        out(j).resp = respj; 
    else
        indrmv = [indrmv j];
    end
 
%             
end

out(indrmv) = [];


%get stim triggered mean and null dist for each site

sitemnboot = cell(1,1);



indk = [];
stim_resp = [];
for k = 1:length(out)
 
       respj = out(k).resp;
       respj = respj - mean(respj(1:baseline));
       stim_resp = [stim_resp; respj];
   
end
sitemn = mean(stim_resp,1);%all the individual stim trials for that site
t = [1:length(sitemn)];
t = t/FS;

%% get the null dist
stimnum = length(out);
Lrsa = length(resp_temp) - indStop_sec*FS;

Nbooti = zeros(stimnum,indStop_sec*FS); %this is the Nboot array of stimnum random samples for each Nboot calc

for k = 1:Nboot
    rr = randi(Lrsa,[1,stimnum]);%get the random sample points for each meand calculation 
    for rri = 1:length(rr)
        respj = resp_temp(rr(rri):rr(rri)+indStop_sec*FS-1);    
        respj = respj - mean(respj(1:baseline));    
        Nbooti(rri,:) = respj;
        
    end

    sitemnboot{1,1} = [sitemnboot{1,1}; mean(Nbooti,1)];%all the boot strapped mean calculation for that site. 
end


mn_boot = mean(sitemnboot{1,1},1);
std_boot = std(sitemnboot{1,1},[],1);


    
    
std_boot_pos = [];
std_boot_neg = [];
p_alpha = 0.025;
    
Nsig = floor((1 - p_alpha)*Nboot);
for bi = 1:length(mn_boot)

    v_bi = sitemnboot{1,1}(:,bi);

    v_bi_ascend = sort(v_bi,'ascend') ;

    v_bi_descend = sort(v_bi,'descend'); 


    std_boot_pos = [std_boot_pos v_bi_ascend(Nsig)];

    std_boot_neg = [std_boot_neg v_bi_descend(Nsig)];
end

    
    
    

            
    


