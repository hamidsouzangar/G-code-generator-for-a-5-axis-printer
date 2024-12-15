% 5-Axis Slicing for FDM with B and C Axis Extensions
% This script simulates 5-axis slicing for a curved beam using a basic
% Cartesian FDM printer with B and C axis extensions.

clc;
clear all;
% Initialization and Input Parameters:
% Clear workspace and figure: Cleans the workspace and creates a new figure.
%Input parameters: 

H_of_layers = 3;        % Layer height
Square_Side = 10;      % Side of the square base
Square_height = 10;    % Height of the square base
Elbow_Degree = 90;     % Angle of the elbow
R_1 = 5;               % Inner radius of the elbow
R_2 = 10;              % Outer radius of the elbow
offset = [1, 0, 0.6];  % Offset of the print bed
nozzle_diameter = 0.5; % Nozzle diameter

% --- Step 1: Square Base Slicing ---
		% Layer Height Calculation for Square Base:
		% Optimal layer height: Calculates the optimal layer height for the square base to minimize excess material.
		% Number of layers: Determines the number of layers for the square base.
	lower_bound = 0.7 * H_of_layers;
	upper_bound = 1.3 * H_of_layers;

	best_number = lower_bound; % Start with the lower bound as the initial best number
	smallest_remainder = mod(Square_Side, best_number);

				for i = lower_bound:0.01:upper_bound;
					remainder = mod(Square_Side, i);
					if remainder < smallest_remainder;
						best_number = i;
						smallest_remainder = remainder;
					end
				end


    H_of_layers_square = best_number;  % layer_height of the first and last cubes

    Layers_Number_square = floor(Square_height/ H_of_layers_square);

% Coordinate System and Visualization Setup:
axis1 = axes('Box', 'on', 'XLim', [-30 50], 'YLim', [-10 50], 'ZLim', [-30 50]);
view(axis1, [45 30]); % Sets the viewing angle.
xlabel('X');	% Adds axis labels and a grid.
ylabel('Y');	% Adds axis labels and a grid.
zlabel('Z');	% Adds axis labels and a grid.
grid on;		
hold(axis1, 'on'); % Keeps the current plot active for subsequent additions.


% Position of the BED coordinate system within the uni-coordinate system
% distance of the first point from orgin of the coordinate
BED_origin = offset; %Defines the origin of the print bed.
BED_origin_2 = BED_origin;  
% Define the vertices of the cube within the BED coordinate system

% pre-initial layer
  vertices =  [
                (BED_origin (1)),                 (BED_origin(2)),                  BED_origin(3);
                (BED_origin (1) - Square_Side),   (BED_origin (2)),                 BED_origin(3);
                (BED_origin (1) - Square_Side),   (BED_origin (2) + Square_Side),   BED_origin(3);
                (BED_origin (1)),                 (BED_origin (2) + Square_Side),   BED_origin(3)
                ];

plot3(vertices(:,1), vertices(:,2), vertices(:,3), 'b');

line([vertices(4,1), vertices(1,1)], [vertices(4,2), vertices(1,2)], [vertices(4,3), vertices(1,3)],  'Color', 'b');

hold on;

% Define the edges of the cube

% Draw the BED coordinate system within the uni-coordinate system
plot3(axis1, [BED_origin(1) BED_origin(1)+50], [BED_origin(2) BED_origin(2)], [BED_origin(3) BED_origin(3)], 'r');
plot3(axis1, [BED_origin(1) BED_origin(1)], [BED_origin(2) BED_origin(2)+50], [BED_origin(3) BED_origin(3)], 'r');
plot3(axis1, [BED_origin(1) BED_origin(1)], [BED_origin(2) BED_origin(2)], [BED_origin(3) BED_origin(3)+50], 'r');

% Draw the process of connecting points to form the cube in BED
% Rotate the BED coordinate system and its contents by 15 degrees around Y-axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%

%% %%%%%%%% Step_two _ elbow %%%%%%%%%%%%%%%%%%%
% Elbow Section Calculations:
Circumference_min = (Elbow_Degree/360) * 2 * R_1 * pi;
Circumference_max = (Elbow_Degree/360) * 2 * (R_2+R_1) * pi;
% calculate H_of_layers for Elbow; % Calculate the 30% range of H_of_layers for elbow
% Similar to the square base, calculates the optimal layer height for the elbow.
	lower_bound = 0.7 * H_of_layers;
	upper_bound = 1.3 * H_of_layers;

	best_number = lower_bound; % Start with the lower bound as the initial best number
	smallest_remainder = mod(Circumference_min, best_number);

			for i = lower_bound:0.01:upper_bound;
				remainder = mod(Circumference_min, i);
				if remainder < smallest_remainder;
					best_number = i;
					smallest_remainder = remainder;
				end
			end


    H_of_layers_Elbow = best_number;  % layer_height of the elbow section;
% Layer Angle and Number of Layers
Gamma_Radian = H_of_layers_Elbow / (R_1 + R_2); %2 * asin ( H_of_layers_Elbow / (2 * (R_1 + R_2)) ) %Gamma = degree of layer %% asin = arcsin
Gamma = Gamma_Radian *(180/pi);



%% Extra Layers and Layer Distribution;
%% No layers; calculating the amounts of bug OR empty eara; it should be almost zero
% Calculates the minimum (incompelete) and maximum (compelete) number of layers for the elbow.
Min_layers = Circumference_min / H_of_layers_Elbow;
Min_layers_number = floor (Min_layers);
Max_layers = Circumference_max / H_of_layers_Elbow;
Max_layers_number = floor (Max_layers);
mode = mod(Max_layers, Min_layers);

All_Layers = Max_layers_number;
final_layers_Number = 0;

if  mode > H_of_layers_Elbow;
    final_layers_Number = floor(mode/H_of_layers_Elbow);
end
% Circumference of sector of a layer = H_of_layers_Elbow
Circumference_final_layers = final_layers_Number * H_of_layers_Elbow;
% nolayers
no_layers = mod(Max_layers, H_of_layers_Elbow);
% degree of no layers area :
nolayers_degree = no_layers / (R_1 + R_2);

All_Layers_Elbow = All_Layers + final_layers_Number;
Gamma_Total = All_Layers_Elbow * Gamma;

Layers_Number_2_Full = floor(Circumference_min/H_of_layers_Elbow); % full layers in elbow
%% calculating final layers of the elbow_shell
Printed_Circumference_min = H_of_layers_Elbow * Layers_Number_2_Full;

Printed_Circumference_max = (Circumference_max * Printed_Circumference_min / Circumference_min);
%%% initial layers = final layers
final_Circumference =  (Circumference_max - Printed_Circumference_max);


%% degree between layers

% rotation degree; it should be almost 90 degree or  Elbow_Degree
rotation_degree = (Elbow_Degree - nolayers_degree) / All_Layers_Elbow;
%% incomplete_ layers numbering;
% number of incomplete spaces,
% spaces between full layers, the printing starts with incomplete layers
incomplete_spaces = Min_layers_number + round((final_layers_Number + 0.1)/(final_layers_Number + 0.3));
% number of incomplete_layers in each incomplete_spaces
incomplete_layers_number = (Max_layers_number - Min_layers_number)/incomplete_spaces;

    All_printing_layers = (2 * Layers_Number_square) + All_Layers_Elbow;
    A = cell(All_printing_layers, 4);

%% Start of the printing simulation %%
%creating matrix that designate which layer must be rotate;
box (1: Layers_Number_square) = zeros; % rotation matrix
box (Layers_Number_square + 1: Layers_Number_square+All_Layers_Elbow) = ones;
box (Layers_Number_square+All_Layers_Elbow+ 1: All_printing_layers ) = zeros;

% box = [0,0,0,1,1,1,1,1,1,1,1]
j  = offset (3); % It must be always zero %% movement of Z in each step
Total_theta = 0; % It must be always zero
counting_1 = 0;
Save_matrix = [0;0;0];
small_theta_beteen = 0 ;   % small_theta_beteen rotated line and actual liNe
%% start the main loop
 for g  = 1:length(box); % number of layers %% length(box)
      A {g,1} = g;
      A {g,3} = 4; % number of line in incompelete layers
      A {g,4} = 0; %to find the rotating layers
  rotation_Degree = box (g)*rotation_degree; % rotation 0 or D-Of_rotation
  Total_theta = Total_theta + rotation_Degree;
%%% if situation for rotation;
% Define the translation to move x=1 axis to origin
translateToOrigin = [1 0 0 -offset(1); 0 1 0 0; 0 0 1 0; 0 0 0 1];

% Define the translation to move back after rotation
translateBack = [1 0 0 offset(1); 0 1 0 0; 0 0 1 0; 0 0 0 1];
if rotation_Degree > 0 % rotation bigger than 0
  cla(axis1); %clear privious axis and its content
  theta = deg2rad(rotation_Degree); % Convert to radians
  rotationMatrixY       = [cos(theta) 0 sin(theta) ;
                           0          1     0      ;
                          -sin(theta) 0 cos(theta)];

% just rotation of BED axis
  theta_2 = deg2rad(Total_theta);
% just rotation of BED axis
  rotationMatrixY_for_BEDAxis = [cos(theta_2) 0 sin(theta_2) ;
                                   0          1     0        ;
                                -sin(theta_2) 0 cos(theta_2)];
% Combine the BED origin translation and rotation into a transformation matrix
origin = [0, 0, 0];
%transformMatrix = [rotationMatrixY, origin'; 0 0 0 1];
% Combine the translations and rotation into a transformation matrix
transformMatrix = translateBack * [rotationMatrixY, [0; 0; 0]; 0 0 0 1] * translateToOrigin;

% just rotation of BED axis
transformMatrix_for_BEDAxis = [rotationMatrixY_for_BEDAxis, BED_origin'; 0 0 0 1];

% Apply the transformation to the BED coordinate system's axes
BED_x_axis_rot = transformMatrix_for_BEDAxis * [40; 0; 0; 1];
BED_y_axis_rot = transformMatrix_for_BEDAxis * [0; 40; 0; 1];
BED_z_axis_rot = transformMatrix_for_BEDAxis * [0; 0; 40; 1];

% Apply the transformation to the vertices of the cube
%rotatedCubeVertices = (transformMatrix * [vertices, ones(size(vertices, 1), 1)]')';





% Position of the BED coordinate system within the uni-coordinate system
BED_origin = offset;

% Draw the rotated BED coordinate system
plot3(axis1, [BED_origin(1) BED_x_axis_rot(1)], [BED_origin(2) BED_x_axis_rot(2)], [BED_origin(3) BED_x_axis_rot(3)], 'g');
plot3(axis1, [BED_origin(1) BED_y_axis_rot(1)], [BED_origin(2) BED_y_axis_rot(2)], [BED_origin(3) BED_y_axis_rot(3)], 'g');
plot3(axis1, [BED_origin(1) BED_z_axis_rot(1)], [BED_origin(2) BED_z_axis_rot(2)], [BED_origin(3) BED_z_axis_rot(3)], 'g');
% Draw the rotated cube in the BED coordinate system



 %this for is for rotating printed layers (past layers)
  for f = 1:g-1 % g = the layer number which is printing;

    X = A{f, 2};
    vertices = X (:,1:3);
    rotatedCubeVertices = (transformMatrix * [vertices, ones(size(vertices, 1), 1)]')';
    vertices = rotatedCubeVertices;
     for i = 1:length(vertices);
        % Plot the current point
        plot3(vertices(i,1), vertices(i,2), vertices(i,3), 'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
        % Connect the point to the previous point if not the first point
         if i > 1;
              line([vertices(i-1,1), vertices(i,1)], [vertices(i-1,2), vertices(i,2)], [vertices(i-1,3), vertices(i,3)], 'Color', 'g', 'LineWidth', 2);
          end

    end
    % Close the square by connecting the last point to the first point

     if A {f,3} == 4;
   line([vertices(end,1), vertices(1,1)], [vertices(end,2), vertices(1,2)], [vertices(end,3), vertices(1,3)], 'Color', 'b', 'LineWidth', 2);
          end

    A{f, 2} = vertices;


  end
end
%%% End of if ;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the equation of a line  (a lyer which in going to be printed fully in the future)
% Calculate the slope (m)
% Ralpha = angle between curren layer to be printed with  the line  (a lyer which in going to be printed fully in the future)

% Elbow_printing
if box (g) == 1;
  counting_1 = counting_1 +1


  %transformation matrix of X and Z after rotation
x_old = Save_matrix (1,1);
z_old = j;

B_rad = - deg2rad( rotation_degree); % Convert to radians

% Define the rotation matrix for the Z-axis
rotation_MatrixY = [cos(B_rad) -sin(B_rad); sin(B_rad) cos(B_rad)];

% Define the old coordinates as a column vector
old_coords = [x_old; z_old];

% Calculate the new coordinates after rotation
new_coords = rotation_MatrixY * old_coords;

% Extract the new x and y coordinates
x_new = new_coords(1);
z_new = new_coords(2);



Save_matrix (:,1) = [new_coords(1);0;new_coords(2)];

j = new_coords(2) + H_of_layers_Elbow;

   % Coordinates for the first line (Line 1)
x1 = 0;
y1 = 0;
x2 = x_new;
y2 = z_new;

% Coordinates for the second line (Line 2)
x3 = 0;
y3 = 0;
x4 = x_new;
y4 = j;

% Direction vectors
A1 = [x2 - x1, y2 - y1];
B1 = [x4 - x3, y4 - y3];

% Dot product
dot_product = dot(A1, B1);

% Magnitudes of the vectors
magnitude_A = norm(A1);
magnitude_B = norm(B1);

% Cosine of the angle
cos_theta = dot_product / (magnitude_A * magnitude_B);

% Angle in radians
theta_rad = acos(cos_theta);

% Convert angle to degrees
small_theta_beteen = abs(rad2deg(theta_rad));


      k = g - Layers_Number_square;
      inn = incomplete_layers_number + 1;
      Ralpha_D = incomplete_layers_number - mod(k-1, inn);

      m = tan(deg2rad( Ralpha_D * rotation_degree)); %slope of the line
      b = -H_of_layers_Elbow;
      z = 0;



% Check if Ralpha is 0
            if Ralpha_D == 0;
                %x = Square_Side
                x = 0;
                X = 0;
                Xt = Square_Side;
            else % Calculate x for other values of Ralpha
            x = (z -b)/m;
            X = x - R_1;
            Xt = Square_Side - X;
            A {g,3} = 3; %to not connect to end points
            end

BED_origin_2 (1) = BED_origin (1) + x_new;

% rotation and printing incomplete layers layers
 vertices =  [
                (BED_origin_2 (1) - Xt),      (BED_origin(2)),                  j;
                (BED_origin_2 (1) ),          (BED_origin (2)),                 j;
                (BED_origin_2 (1) ),          (BED_origin (2) + Square_Side),   j;
                (BED_origin_2 (1) - Xt),      (BED_origin (2) + Square_Side),   j
                ];
%rotatedCubeVertices = (transformMatrix * [vertices, ones(size(vertices, 1), 1)]')';
 %          vertices = rotatedCubeVertices

         A {g,2} = vertices;
     for i = 1:length(vertices);
        % Plot the current point
        plot3(vertices(i,1), vertices(i,2), vertices(i,3), 'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'g')
        pause(0.1); % Pause for 1 second to visualize the process

        % Connect the point to the previous point if not the first point
        if i > 1
            line([vertices(i-1,1), vertices(i,1)], [vertices(i-1,2), vertices(i,2)], [vertices(i-1,3), vertices(i,3)], 'Color', 'r', 'LineWidth', 2);
        end
        pause(0.1);
     end
        % Close the square by connecting the last point to the first point
        if Ralpha_D == 0;
         line([vertices(end,1), vertices(1,1)], [vertices(end,2), vertices(1,2)], [vertices(end,3), vertices(1,3)], 'Color', 'b', 'LineWidth', 2);
       end

      pause(0.1)
         else
     j = j + H_of_layers_square;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% printing layers_square
if box (g) == 0;
  counting_22 = g
  vertices =  [
              (BED_origin_2 (1)),                (BED_origin(2)),                  j;
              (BED_origin_2 (1) - Square_Side),  (BED_origin (2)),                 j;
              (BED_origin_2 (1) - Square_Side),  (BED_origin (2) + Square_Side),   j;
              (BED_origin_2 (1)),                (BED_origin (2) + Square_Side),   j
                ];
%rotatedCubeVertices = (transformMatrix * [vertices, ones(size(vertices, 1), 1)]')';
 %          vertices = rotatedCubeVertices
     for i = 1:length(vertices);
        % Plot the current point
        plot3(vertices(i,1), vertices(i,2), vertices(i,3), 'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
        pause(0.1); % Pause for 1 second to visualize the process

        % Connect the point to the previous point if not the first point
        if i > 1
            line([vertices(i-1,1), vertices(i,1)], [vertices(i-1,2), vertices(i,2)], [vertices(i-1,3), vertices(i,3)], 'Color', 'b', 'LineWidth', 2);
        end
     end
        % Close the square by connecting the last point to the first point
         line([vertices(end,1), vertices(1,1)], [vertices(end,2), vertices(1,2)], [vertices(end,3), vertices(1,3)], 'Color', 'b', 'LineWidth', 2);
         A {g,2} = vertices;
         pause(0.1)
end
end




view(axis1, [45 30]);

hold(axis1, 'off');
