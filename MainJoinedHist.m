function [] = MainJoinedHist( y1, y2, y3, Text, fileName, IsWorstShown)
z1 = sort(y1);
z2 = sort(y2);
z3 = sort(y3);

n = length(z1);
%     x1 = z1;
%     x2 = z2;
%     x3 = z3;
 x1 = z1(round(0.02*n)+1:round(0.98*n),1);
 x2 = z2(round(0.02*n)+1:round(0.98*n),1);
 x3 = z3(round(0.02*n)+1:round(0.98*n),1);
 
figure1 = figure;

set(figure1,'defaulttextinterpreter','latex');
set(gca, 'FontSize', 25);

[f1, xi1] = ksdensity(x1);
[f2, xi2] = ksdensity(x2);
[f3, xi3] = ksdensity(x3);

graph1 = plot(xi1, f1, 'Color', [0, 0, .7]);
hold on
graph2 = plot(xi2, f2, 'Color', [0, .7, 0]);



legend ('Normal', 'Extreme')

if (IsWorstShown == 1)
    hold on
    graph3 = plot(xi3, f3, 'Color', [.7, 0, 0]);
    medianLine (x3, xi3, f3,[.7, 0, 0])
    set (graph3, 'LineWidth', 2);
    legend ('Normal', 'Extreme', '10Extreme')
end

h_legend= legend ('boxoff');
set(h_legend,'FontSize', 24);

 
 set (graph1, 'LineWidth', 2);
 set (graph2, 'LineWidth', 2);

 medianLine (x1, xi1, f1, [0, 0, .7])
 medianLine (x2, xi2, f2,[0, .7, 0])

 title( Text,'Interpreter','latex', 'fontsize', 20);
h = get(gca, 'title');
hold off

set(h_legend,'FontSize', 24);
set(gcf, 'Position', get(0,'Screensize'));
set(gcf, 'PaperPositionMode', 'auto')
saveas(gcf, fileName, 'png')