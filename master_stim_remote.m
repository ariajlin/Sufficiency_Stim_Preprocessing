function master_stim_remote
%To run on the server: (1) Upload Trial Notes to data_store1, (2) Hilbert transform
%relevant blocks, (3) Upload stim parameter list to userdata

%Fill in:
pN = ; %Number only
blocks = []; %Number array
goodIndices = []; %Number array
yMin = ; %Use to specify plotting window
yMax = ; %Use to specify plotting window

useremotionpath = getenv('EMOTION_ROOT');
addpath(genpath(fullfile(useremotionpath,'Code','Emotion_Preprocessing')));
PatientID = ['EC' num2str(pN)];
getSufficiencyStimStartTimes_remote(pN, blocks, goodIndices); %I believe this will save output files to userdata directory. In that case, you'll need to manually move them to data_store1 in a "SufficiencyStim" folder (see other pts for formatting) before running the remaining functions.

%After completing and checking output files make sense, run the following:
processingInstructions_stim_mod(pN, blocks); %This function call requires a patient-specific script, only inputs you need to change are Patient ID and blocks

%To plot EDA analysis:
stimERPs_withBootstrap_v4_mod(PatientID,'phasic',yMin,yMax,15); %15s plotting window
stimERPs_withBootstrap_v4_mod(PatientID,'phasic',yMin,yMax,17); %17s plotting window
stimERPs_withBootstrap_v4_mod(PatientID,'phasic',yMin,yMax,20); %20s plotting window

%To plot RSA analysis:
% stimERPs_withBootstrap_v4_mod(PatientID,'RSA',yMin,yMax,15,8);
% stimERPs_withBootstrap_v4_mod(PatientID,'RSA',yMin,yMax,15,16);
% stimERPs_withBootstrap_v4_mod(PatientID,'RSA',yMin,yMax,17,8);
% stimERPs_withBootstrap_v4_mod(PatientID,'RSA',yMin,yMax,17,16);
% stimERPs_withBootstrap_v4_mod(PatientID,'RSA',yMin,yMax,20,8);
% stimERPs_withBootstrap_v4_mod(PatientID,'RSA',yMin,yMax,20,16);


