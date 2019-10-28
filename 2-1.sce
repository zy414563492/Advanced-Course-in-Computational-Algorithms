figure(1);
clf;
//2 6
a = 5;
b = 1.5;

t = linspace(0, 12*%pi, 300);

//x = (a - b) * cos(t) + b * cos((a - b) / b * t);
//y = (a - b) * sin(t) - b * sin((a - b) / b * t);

x = (a - b) * cos(3 * t) + b * cos((a - b) / b * 3 * t);
y = (a - b) * sin(2 * t) - b * sin((a - b) / b * 2 * t);

plot(x, y);


