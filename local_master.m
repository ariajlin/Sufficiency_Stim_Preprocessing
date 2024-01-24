pNs = {'EC280','EC281','EC286','EC288','EC292','EC293','EC296'};

%For plotting EDA:
yMins = [-2.5,-2,-3,-0.5,-1,-2,-0.3];
yMaxs = [4,4.5,3,1,6,2,0.3];

for i = 1:length(pNs)
    p = string(pNs(i))
    stimERPs_withBootstrap_v4(p,'phasic',yMins(i),yMaxs(i),15); %15s plot window
    stimERPs_withBootstrap_v4(p,'phasic',yMins(i),yMaxs(i),17); %17s plot window
    stimERPs_withBootstrap_v4(p,'phasic',yMins(i),yMaxs(i),20); %20s plot window
end

%For plotting RSA:
yMins = [-3,-2,-2,-2,-2.5,-1,-2];
yMaxs = [4,5,3,2,4,6,2];

for i = 1:length(pNs)
    stimERPs_withBootstrap_v4(pNs(i),'phasic',yMins(i),yMaxs(i),15,8); %15s plot window, 8s melSpec window
    pNs(i)
    stimERPs_withBootstrap_v4(pNs(i),'phasic',yMins(i),yMaxs(i),15,16); %15s plot window, 16s melSpec window
    stimERPs_withBootstrap_v4(pNs(i),'phasic',yMins(i),yMaxs(i),17,8); %17s plot window, 8s melSpec window
    stimERPs_withBootstrap_v4(pNs(i),'phasic',yMins(i),yMaxs(i),17,16); %17s plot window, 16s melSpec window
    stimERPs_withBootstrap_v4(pNs(i),'phasic',yMins(i),yMaxs(i),20,8); %20s plot window, 8s melSpec window
    stimERPs_withBootstrap_v4(pNs(i),'phasic',yMins(i),yMaxs(i),20,16); %20s plot window, 16s melSpec window
end
