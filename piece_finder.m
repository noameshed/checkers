%% Image Processing Final Project
% Noam Eshed
% 4/29/2018

close all; clc; clear all;
% Path to image
im = (imread('C:\Users\noam_\Documents\RPI\Image Processing\Final Project\Checkers\top4.jpg'));

%% Part 1: Transform the board into a square
% Transformation
% Find the corners
% Click on corners in order: top left, top right, bottom left, bottom right
imshow(im)
[xi, yi] = getpts;
close all

u1 = [xi(1) yi(1)];
u2 = [xi(2) yi(2)];
u3 = [xi(3) yi(3)];
u4 = [xi(4) yi(4)];

% Moving Points
U = [u1; u2; u3; u4];

% Fixed Points
side = 1000;
X = [0 0;...
    side 0;...
    0 side;...
    side side];

% Make transformation object to project corners from U to X
TF = fitgeotrans(U,X,'projective') ;
[Square, obj] = imwarp(im, TF);
%figure; imshow(Square)

% Crop transformed image to shape
Sq1 = Square(abs(round(obj.YWorldLimits(1)))+[0:side],...
    abs(round(obj.XWorldLimits(1)))+[0:side],:);
figure; imshow(Sq1)

%% Part 2: Thresholding and finding Canny edges
% Method 1 - RGB2GRAY
gs1 = rgb2gray(Sq1);

% Method 2 - Red Channel
gs2 = Sq1(:,:,1); 

% Get canny edges
BW1 = edge(gs1, 'canny');

T = graythresh(gs2);
%BW2 = edge(gs2, 'canny');   % has more complete outside borders
BW2 = edge(gs2<255*T, 'Canny');
%imshow(BW2-BW1);

figure;
set(gcf, 'Position', get(0, 'Screensize'));
subplot(1,4,1), imshow(gs2); title('red channel')
subplot(1,4,2), imshow(BW2); title('Canny Edges')

% Morphological opening
h = strel('line', 4, 0);
v = strel('line', 4, 90);
o1 = imopen(BW2, h);
o2 = imopen(BW2, v);
o = o1+o2;

% Opening a second time w/ a larger strel
h = strel('line', 6, 0);
v = strel('line', 6, 90);
o3 = imopen(o, h);
o4 = imopen(o, v);
fin_o = o3 + o4;

subplot(1,4,3), imshow(o); title('With Morph Opening')
subplot(1,4,4), imshow(fin_o); title('Second round opening')

%% Finding Hough edges
[H, theta, rho] = hough(fin_o);
peaks = houghpeaks(H, 30);
lines = houghlines(fin_o, theta, rho, peaks);

% Plot the Hough lines
figure; imshow(Sq1); hold on
for i = 1:length(lines)
   xy = [lines(i).point1; lines(i).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
end
hold off
%% Combine Hough line segments which correspond to the same line
x = zeros(2,length(lines));     % all vertical segments
mx = zeros(1,length(lines));    % their slopes
y =  zeros(2,length(lines));    % all horizontal segments
my =  zeros(1,length(lines));   % their slopes

for i = 1:length(lines)
    x1 = lines(i).point1(1);
    y1 = lines(i).point1(2);
    x2 = lines(i).point2(1);
    y2 = lines(i).point2(2);
    
    m = (y2-y1)/(x2-x1);            % slope of the line   
 
    if (abs(m) > 60) || (m == Inf)       % Vertical Line
        % Check if a similar segmant has already been included
        for j = 1:i
            % For all segments already included
            if (i ~= j) && (abs(mean(x1, x2) - mean(x(:,j))) < 15) %&& ~isequal(x(:,j), [0;0])
                % If two segments are too similar, don't include the new
                % one
                x(:,i) = [0,0]; break
            else
                % Add new segment if it's far enough away
                x(:,i) = [x1; x2];
                mx(i) = m;               
            end
        end        
        
    elseif abs(m) < 25       % Horizontal Line
        for j = 1:i
            % For all segments already included
            if (i ~= j) && (abs(mean(y1, y2) - mean(y(:,j))) < 15) && ~isequal(y(:,j),[0;0])
                % If two segments are too similar, don't include the new
                % one
                y(:,i) = [0;0]; break
            else
                % Add new segment if it's far enough away
                y(:,i) = [y1; y2];
                my(i) = m;               
            end
        end       
    end    
end

%% First remove any points which aren't actually lines ([0;0])
xprime = [];
for z = 1:length(x) % Plot vertical lines
    if x(:,z) > [0;0]
        xprime = [xprime, x(:,z)];
    end
end
yprime = [];
for z = 1:length(y) % Plot horizontal lines
    if y(:,z) > [0;0]
        yprime = [yprime, y(:,z)];
    end
end

xprime = sort(xprime, 2);
yprime = sort(yprime, 2);

%% Check that image edges are in and add them if not
% If we didn't detect any of the outer edges, add them in:
if size(xprime,2) <2 || mean(xprime(1)) > 20
    xprime = [[1;1], xprime];
end
if size(xprime,2) <2 || mean(xprime(:,(size(xprime, 2)))) <980
    xprime = [xprime, [size(Sq1, 1);size(Sq1, 1)],];
end
if size(yprime,2) <2 || mean(yprime(1)) > 20
    yprime = [[1;1], yprime];
end
if size(yprime,2) <2 || mean(yprime(:,(size(yprime, 2)))) <980
    yprime = [yprime, [size(Sq1, 1);size(Sq1, 1)],];
end

%% If we don't have 9 lines in a certain orientation, interpolate the
% Note that this won't work if Hough lines doesn't pick up enough
% edges. This supposes that most of the edges are found (i.e. the
% "mode" of distances is the actual distance between two lines on the
% board)
figure; imshow(Sq1); hold on
while size(xprime,2) < 9
    distx = [xprime [0;0]]-[[0;0] xprime];
    distx = distx(:,2:size(distx,2)-1);
    round_distx = round(distx/15)*15;    
    modex = mode(round_distx(:)); 
    for i = 1:size(round_distx,2)
        if abs(mean(distx(:,i)) - modex) > 25
            
            % If we're missing lines somewhere
            prevline = xprime(:,i);
            nextline = xprime(:,i+1);
            top_diff = prevline(1)-nextline(1);
            bottom_diff = prevline(2)-nextline(2);
            
            % How many lines do we need to add?
            num_lines = round(abs(top_diff/modex));
            
            % Add the line(s)
            for j = 1:num_lines-1
                top_diff = nextline(1)-prevline(1);
                bottom_diff = nextline(2)-prevline(2);
                num_lines = round(abs(top_diff/modex));
                newline = prevline + [round(top_diff/num_lines);round(bottom_diff/num_lines)];
                
                prevline = newline;
                nextline = xprime(:,i+j); 
                
                plot(newline, [0;size(Sq1,1)],'LineWidth',2,'Color','blue');
                xprime = [xprime(:, 1:i+j-1),newline,xprime(:,i+j:length(xprime))];
            end
        end
    end
end
while size(yprime,2) < 9
  
    disty = [yprime [0;0]]-[[0;0] yprime];
    disty = disty(:,2:size(disty,2)-1);
    round_disty = round(disty/15)*15;    
    modey = mode(round_disty(:));
    
    for i = 1:size(round_disty,2)
        if abs(mean(disty(:,i)) - modey) > 25
            
            % If we're missing lines somewhere
            prevline = yprime(:,i);
            nextline = yprime(:,i+1);
            top_diff = prevline(1)-nextline(1);
            bottom_diff = prevline(2)-nextline(2);
            
            % How many lines do we need to add?
            num_lines = round(abs(top_diff/modey));
            
            % Add the line(s)
            
            for j = 1:num_lines-1
                top_diff = nextline(1)-prevline(1);
                bottom_diff = nextline(2)-prevline(2);
                num_lines = round(abs(top_diff/modey));
                newline = prevline + [round(top_diff/num_lines);round(bottom_diff/num_lines)];
                
                prevline = newline;
                nextline = yprime(:,i+j); 
                
                plot([0;size(Sq1,1)],newline, 'LineWidth',2,'Color','blue');
                yprime = [yprime(:, 1:i+j-1),newline,yprime(:,i+j:length(yprime))];
            end
        end
    end
end
hold off
%% Plot the interpolated lines
figure; imshow(Sq1); hold on
for z = 1:length(xprime) % Plot vertical lines
    plot(xprime(:,z), [0;size(Sq1,1)],'LineWidth',2,'Color','green');
end

for z = 1:length(yprime) % Plot horizontal lines
    plot([0;size(Sq1,2)], yprime(:,z),'LineWidth',2,'Color','green');
end

%% Thresholding Methods
% Original (red channel)
%figure
%subplot(2,3,1), imshow(Sq1), title('Original')

% 1. Red Channel
gs = Sq1(:,:,1); 
%subplot(2,3,2); imshow(gs), title('Grayscale')

% 2. Global Otsu thresholding
T = graythresh(gs);
global_otsu = gs<255*T;
%subplot(2,3,3); imshow(global_otsu, []), title('Global thresholding-Otsu')

% 3. Regional thresholding using Otsu
fun = @(block_struct) otsu_local_thresh(block_struct.data);
B = blockproc(gs, [150 150], fun);
%subplot(2,3,4); imshow(B, []), title('Regional Thresholding - Otsu')

% 4. Simple thresholding
binary = (gs)<55;   % FUTURE WORK - Play around with this number
%subplot(2,3,5);  imshow(binary), title('Simple thresholding')

% 5. Find Canny edges
BW = edge(gs<255*T, 'Canny');   % can add threshold value
%subplot(2,3,6); imshow(BW), title('Canny Edges')

radius = [40 60];
% look for bright circles
T = 0;
[centers1, radii1] = imfindcircles(Sq1, radius, 'EdgeThreshold', T);
[centers2, radii2] = imfindcircles(gs, radius, 'EdgeThreshold', T);
[centers3, radii3] = imfindcircles(global_otsu, radius, 'EdgeThreshold', T);
[centers4, radii4] = imfindcircles(B, radius, 'EdgeThreshold', T);
[centers5, radii5] = imfindcircles(binary, radius, 'EdgeThreshold', T);
[centers6, radii6] = imfindcircles(BW, radius, 'EdgeThreshold', T);

% look for dark circles
[centers, radii] = imfindcircles(Sq1, radius, 'EdgeThreshold', T, 'ObjectPolarity', 'dark');
centers1 = [centers1; centers]; radii1 = [radii1; radii];
[centers, radii] = imfindcircles(gs, radius, 'EdgeThreshold', T, 'ObjectPolarity', 'dark');
centers2 = [centers2; centers]; radii2 = [radii2; radii];
[centers, radii] = imfindcircles(global_otsu, radius, 'EdgeThreshold', T, 'ObjectPolarity', 'dark');
centers3 = [centers3; centers]; radii3 = [radii3; radii];
[centers, radii] = imfindcircles(B, radius, 'EdgeThreshold', T, 'ObjectPolarity', 'dark');
centers4 = [centers4; centers]; radii4 = [radii4; radii];
[centers, radii] = imfindcircles(binary, radius, 'EdgeThreshold', T, 'ObjectPolarity', 'dark');
centers5 = [centers5; centers]; radii5 = [radii5; radii];
[centers, radii] = imfindcircles(BW, radius, 'EdgeThreshold', T, 'ObjectPolarity', 'dark');
centers6 = [centers6; centers]; radii6 = [radii6; radii];


figure;
subplot(2,3,1); imshow(Sq1); viscircles(centers1, radii1,'EdgeColor','b');title('Original')
subplot(2,3,2); imshow(gs); viscircles(centers2, radii2,'EdgeColor','b');title('Grayscale')
subplot(2,3,3); imshow(global_otsu); viscircles(centers3, radii3,'EdgeColor','b');title('Global thresholding-Otsu')
subplot(2,3,4); imshow(B); viscircles(centers4, radii4,'EdgeColor','b');title('Regional Thresholding - Otsu')
subplot(2,3,5); imshow(binary); viscircles(centers5, radii5,'EdgeColor','b'); title('Simple thresholding')
subplot(2,3,6); imshow(BW); viscircles(centers6, radii6,'EdgeColor','b');title('Canny Edges')


%% Store circle locations and radii

xycenters_raw = [centers1;centers2;centers3;centers4;centers5;centers6];
xyradii_raw = [radii1;radii2;radii3;radii4;radii5;radii6];

%% Add the circles to the original image combining these methods

[val, idx] = max([size(centers1,1), size(centers2,1), size(centers3,1), size(centers4,1), size(centers5,1), size(centers6,1)]);
figure; imshow(Sq1), hold on
if idx == 1 &&  ~isempty(centers1)
    viscircles(centers1, radii1,'EdgeColor','b');
    %xc = centers1(:,1); yc = centers1(:,2);
elseif idx == 2 &&  ~isempty(centers2)
    viscircles(centers2, radii2,'EdgeColor','b');
    %xc = centers2(:,1); yc = centers2(:,2);
elseif idx == 3 &&  ~isempty(centers3)
    viscircles(centers3, radii3,'EdgeColor','b');
    %xc = centers3(:,1); yc = centers3(:,2);
elseif idx == 4 &&  ~isempty(centers4)
    viscircles(centers4, radii4,'EdgeColor','b');
    %xc = centers4(:,1); yc = centers4(:,2);
elseif idx == 5 &&  ~isempty(centers5)
    viscircles(centers5, radii5,'EdgeColor','b');
    %xc = centers5(:,1); yc = centers5(:,2);
elseif idx == 6 &&  ~isempty(centers6)
    viscircles(centers6, radii6,'EdgeColor','b');
    %xc = centers6(:,1); yc = centers6(:,2);
end

%% Combine circles found by these methods, remove/average
% circles that are too close
xycenters_raw = sortrows(xycenters_raw);
xycenters = [];
xyradii = [];
j = 2;
for i = 1:size(xycenters_raw,1)-1
    add = xycenters_raw(i,:);
    addr = xyradii_raw(i,:);
    next_pt = xycenters_raw(j,:);
    next_r = xyradii_raw(j,:);

    if  abs(mean(add - next_pt)) < 10
        % If we already have a close point, just average the coordinatess
        j = j+1;
        add = [(add(1)+next_pt(1))/2.0, (add(2)+next_pt(2))/2.0];
        addr = (addr + next_r)/2;
    else
        % If this is a new point, then push i up to our current index
        j = i;
        xycenters = [xycenters; add];
        xyradii = [xyradii; addr];
    end
    
end

%figure; imshow(Sq1); hold on

figure; imshow(Sq1); hold on
for z = 1:length(xprime) % Plot vertical lines
    plot(xprime(:,z), [0;size(Sq1,1)],'LineWidth',2,'Color','green');
end
for z = 1:length(yprime) % Plot horizontal lines
    plot([0;size(Sq1,2)], yprime(:,z),'LineWidth',2,'Color','green');
end
viscircles(xycenters, xyradii, 'EdgeColor', 'b');


%% Find which squares have pieces

% The 'name' matrix is simply for returning the name of the piece location 
name = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8';...
               'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8';...
               'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8';...
               'D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8';...
               'E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8';...
               'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8';...
               'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8';...
               'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8'];
col = 1:8;
row = ['A','B','C','D','E','F','G','H'];
% all_pcs is an 8x8 matrix that gives the locations of all pieces
% level 1 is x coordinate, level 2 is y coordinate
% if the values are 0, there is no piece there
all_pcs = zeros(8,8,2);         
all_locs = [];

% Loop through all of the pieces and find which square it is in
for i = 1:size(xycenters,1)

    % for center point i, find its corresponding row and column
    xc = xycenters(i,1);
    yc = xycenters(i,2);
    xloc = 0;
    
    for linex = 1:size(xprime,2)
        if xc > mean(xprime(:,linex)) && xc < mean(xprime(:,linex+1))
            xloc = linex;
        end
    end
    
    yloc = 0;
    for liney = 1:size(yprime,2)
        if (yc > mean(yprime(:,liney))) && (yc < mean(yprime(:,liney+1)))            
            yloc = liney;
        end
    end
    %xloc, yloc
    all_pcs(yloc, xloc, :) = [xc, yc];
    all_locs = [all_locs;(strcat(row(yloc), num2str(col(xloc))))];
end
disp('There are checker pieces at locations:')
disp(sortrows(all_locs))

