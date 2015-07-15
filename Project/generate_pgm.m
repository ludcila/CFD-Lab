clear all; clc; close all;

[folder, ~, ~] = fileparts(which('generate_pgm'))

C_B = 0;
C_F = 1;

% Box

imax = 100;
jmax = 20;

domain = ones(jmax, imax);


% % Square obstacle
% domain(9:10, 12) = C_B;
% domain(9:11, 11) = C_B;
% domain(10:12, 10) = C_B;
% domain(11:12, 9) = C_B;
% imshow(domain);
% imwrite(domain, [folder, '/karman_vortex_street.pgm'], 'encoding', 'ASCII', 'maxvalue', 1);

% % Flow over step
% domain(jmax/2+1:end, 1:jmax/2) = C_B;
% imshow(domain);
% imwrite(domain, [folder, '/flow_over_step.pgm'], 'encoding', 'ASCII', 'maxvalue', 1);

% % Flow over step
% imwrite(domain, [folder, '/plane_shear_flow.pgm'], 'encoding', 'ASCII', 'maxvalue', 1);

% Driven cavity
imax = 100;
jmax = 40;
[X Y] = meshgrid(linspace(0, 100), linspace(0, 40, 40));
domain = zeros(jmax, imax);

% center dam
% domain(2:end, 2*imax/5+1:3*imax/5) = 1;

% drop
% % domain(5:15, 45:55) = 1;
% domain = domain + sqrt((X-50).^2 + (Y-10).^2) < 5;
% domain(30:end, :) = 1;

% bubble
% domain(10:end, :) = 1;
% % domain(25:35, 45:55) = 0;
% domain = domain - (sqrt((X-50).^2 + (Y-20).^2) < 1);

% Dam break
imax = 100;
jmax = 40;
% domain(2:end, 1:imax/5) = 1;
domain(2:end, 2*imax/5+1:3*imax/5) = 1;
imwrite(domain, ['dam_break_center.pgm'], 'encoding', 'ASCII', 'maxvalue', 1);

% Dam break on step
imax = 50;
jmax = 50;
domain = zeros(jmax, imax);
domain(31:50, 1:20) = 2;
domain(10:30, 1:20) = 1;
% imwrite(domain, ['dam_break_step.pgm'], 'encoding', 'ASCII', 'maxvalue', 1);

% Dam break small obstacle
imax = 50;
jmax = 50;
domain = zeros(jmax, imax);
domain(47:50, 30:33) = 2;
domain(11:end, 1:15) = 1;
% imwrite(domain, ['dam_break_obstacle.pgm'], 'encoding', 'ASCII', 'maxvalue', 1);


% domain(20:end, :) = 1;

% domain(5,5) = 1;
