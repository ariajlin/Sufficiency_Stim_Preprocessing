function p = stdshade_MeanSem_v2(amean,astd,alpha,acolor,F,smth)
% usage: stdshading(amatrix,alpha,acolor,F,smth)
% plot mean and sem/std coming from a matrix of data, at which each row is an
% observation. sem/std is shown as shading.
% - acolor defines the used color (default is red) 
% - F assignes the used x axis (default is steps of 1).
% - alpha defines transparency of the shading (default is no shading and black mean line)
% - smth defines the smoothing factor (default is no smooth)
% smusall 2010/4/23

if exist('acolor','var')==0 || isempty(acolor)
    acolor='r'; 
end

if(size(amean,1)>size(amean,2))
    amean = amean';
    astd = astd';
end

if exist('F','var')==0 || isempty(F); 
    F=1:max(size(amean));
end

if exist('smth','var'); if isempty(smth); smth=1; end
else smth=1;
end  

if ne(size(F,1),1)
    F=F';
end
% 
% amean=smooth(nanmean(amatrix),smth)';
% %astd=nanstd(amatrix); % to get std shading
% astd=nanstd(amatrix)/sqrt(size(amatrix,1)); % to get sem shading
if length(astd) > 1
    if exist('alpha','var')==0 || isempty(alpha) 
        fill([F fliplr(F)],[amean+2*astd fliplr(amean-2*astd)],acolor,'linestyle','none');
        acolor='k';
    else
        fill([F fliplr(F)],[amean+2*astd fliplr(amean-2*astd)],acolor, 'FaceAlpha', alpha,'linestyle','none');    
    end
end

if ishold==0
    check=true; else check=false;
end

hold on;
p = plot(F,amean,'Color',acolor,'linewidth',1.5); %% change color or linewidth to adjust mean line

if check
    hold off;
end

end



