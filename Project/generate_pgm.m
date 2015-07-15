clear all; clc; close all;

% =========================================================================

% Options:

CENTER_DAM = 1;
LEFT_DAM = 2;
DAM_STEP = 3;
DAM_OBSTACLE = 4;
DROP = 5;
BUBBLE = 5;

% =========================================================================

% Specify scenario to generate:

scenario = DAM_OBSTACLE;

% =========================================================================




if scenario == CENTER_DAM
    
    imax = 125;
    jmax = 50;
    domain = zeros(jmax, imax);
    domain(2:end, 2*imax/5+1:3*imax/5) = 1;
    name = 'dam_break_center.pgm';
    
elseif scenario == LEFT_DAM
    
    imax = 125;
    jmax = 50;
    domain = zeros(jmax, imax);
    domain(2:end, 1:imax/5) = 1;
    name = 'dam_break_left.pgm';
    
elseif scenario == DAM_STEP
    
    imax = 125;
    jmax = 50;
    domain = zeros(jmax, imax);
    domain(jmax/2+1:jmax, 1:imax*2/5) = 2;
    domain(13:jmax/2, 1:imax*2/5) = 1;
    name = 'dam_break_step.pgm';
    
elseif scenario == DAM_OBSTACLE
    
    imax = 80;
    jmax = 80;
    domain = zeros(jmax, imax);
    domain(75:80, 39:41) = 2;
    domain(41:end, 1:20) = 1;
    name = 'dam_break_obstacle.pgm';
    
elseif scenario == DROP
    
    imax = 125;
    jmax = 50;
    domain = zeros(jmax, imax);
    [X Y] = meshgrid(linspace(1, imax, imax), linspace(1, jmax, jmax));
    domain = domain + sqrt((X-imax/2).^2 + (Y-12.5).^2) < 3;
    domain(38:end, :) = 1;
    name = 'drop.pgm';
    
elseif scenario == BUBBLE
    domain(10:end, :) = 1;
    % domain(25:35, 45:55) = 0;
    domain = domain - (sqrt((X-50).^2 + (Y-20).^2) < 1);
    name = 'bubble.pgm';
end

imshow(domain)
[folder, ~, ~] = fileparts(which('generate_pgm'));
imwrite(domain, [folder,'/',name], 'encoding', 'ASCII', 'maxvalue', 1);

