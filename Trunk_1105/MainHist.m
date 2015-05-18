function [] = MainHist( y, nbins, Text, fileName)
z = sort(y);
n = length(z);
%    x = z;
 x = z(round(0.02*n)+1:round(0.98*n),1);
figure1 = figure;

% figure('units','normalized','outerposition',[0 0 1 1])
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 17.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])


set(figure1,'defaulttextinterpreter','latex');
%histfit(x, nbins, 'kernel');
set(gca, 'FontSize', 25);
[heights,centers] = hist(x, nbins, 'kernel', 'normal');
heights = heights/trapz(centers, heights);

bar(centers,heights, 'FaceColor',[0.7,0.7, 0.7], 'EdgeColor',[0.7,0.7, 0.7]);
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
try
fill(t([p p p+1 p+1]),[0 heights([p p]),0], 'r')
catch
end

i = 1;
while (centers(i) < quantile(x, 0.50))&&(i < length(centers))
    i=i+1;
end
if centers(i) - 0 < w/2
    p = i;
else
    p = i-1;
end
fill(t([p p p+1 p+1]),[0 heights([p p]),0], [0, 0, .7])

[f, xi] = ksdensity(x);
graph = plot(xi, f, 'Color', [0, 0, .7]);
set (graph(1), 'LineWidth', 4.5);
title( Text,'Interpreter','latex', 'fontsize', 40);
h = get(gca, 'title');
set(h, 'FontName', 'Helvetica', 'fontsize', 40, 'FontWeight','bold')
hold off

set(gcf, 'Position', get(0,'Screensize'));
set(gcf, 'PaperPositionMode', 'auto')
saveas(gcf, fileName, 'jpg')