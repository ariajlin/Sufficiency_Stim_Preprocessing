function stdshade_mn_sem_computed_shadeOnly(x,ymean,ysem_pos,ysem_neg, rgb,alpha,lineW)

%modified to plot the shade only, no line

%plots the shaded the mean and sem when the ymean and ysem have already been calculated
%i.e. it shades the area of ymean + ysem to ymean - sem with color rgb and
%transparencey of alpha. 
%lineW = line width. 
%alpha is the transparence of the shading (0-1). 
%PH 2021


hold on
curve1 = ymean + ysem_pos;
curve2 = ymean + ysem_neg;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
hfill = fill(x2, inBetween,'g');
hfill.FaceColor = rgb;
hfill.EdgeColor = rgb;
hfill.FaceAlpha = alpha;
hfill.EdgeAlpha = alpha;
%plot(x, ymean,'Color',rgb, 'LineWidth', lineW);
hold off
