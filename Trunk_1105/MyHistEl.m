function [] = MyHistEl( y, target, nbins, Text)
%MY. Summary of this function goes here
%  Detailed explanation goes here

if target ~= 0
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
else
    z = sort(y);
    n = length(z);
    x = z;
%     x = z(round(0.01*n)+1:round(0.99*n),1);
    figure1 = figure;
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
    fill(t([p p p+1 p+1]),[0 heights([p p]),0], 'r')
    
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
end

end

