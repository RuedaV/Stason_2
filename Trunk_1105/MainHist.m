function [] = MainHist( y, nbins, Text, fileName)

x = (y-target)/target;
figure1 = figure;
set(figure1,'defaulttextinterpreter','latex');
%histfit(x, nbins, 'kernel');
set(gca, 'FontSize', 24);
[heights,centers] = hist(x, nbins, 'kernel', 'normal');
heights = heights/trapz(centers, heights);
bar(centers,heights);
hold on
n = length(centers);
w = centers(2)-centers(1);
t = linspace(centers(1)-w/2,centers(end)+w/2,n+1);
i =1;
while (centers(i) < 0)&&(i < length(centers))
    i=i+1;
end
if centers(i) - 0 < w/2
    p = i;
else
    p = i-1;
end
fill(t([p p p+1 p+1]),[0 heights([p p]),0], [1, 1, 0])
[f, xi] = ksdensity(x);
graph = plot(xi, f, 'r');
set (graph(1), 'LineWidth', 2.5);
title( Text,'Interpreter','latex', 'fontsize',24);
hold off

saveas(gcf, fileName, 'png')