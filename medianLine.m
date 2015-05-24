function [] = medianLine (x1, xi1, f1, color)

m1 = median(x1);
n = length(xi1);
w = xi1(2)-xi1(1);
i = 1;
while (xi1(i) < m1)&&(i < n)
    i=i+1;
end
if xi1(i) - m1 < w/2
    p = i;
else
    p = i-1;
end
x = [m1, m1];
y = [0, f1(p)];
line(x, y, 'Color', color, 'LineStyle', '--', 'LineWidth', 1.5);