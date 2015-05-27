clear all; close all;

N = 50;
k = 0:N-1;
h = 1/(N-1);
ws = linspace(0, 2, 100);

err_reds = zeros(1, 10);
figure; hold on;
i = 1;
for w = ws
    err_red = (1-w) + w*cos(pi*k*h);
    plot(k, err_red);
    err_reds(i) = sum(err_red.^2);
    i = i + 1;
end
ylim([-3, 3]);
grid on;
figure();
plot(ws, sqrt(err_reds)/20);
title('Weighted Jacobi');
ylabel('Error reduction');
xlabel('\omega');