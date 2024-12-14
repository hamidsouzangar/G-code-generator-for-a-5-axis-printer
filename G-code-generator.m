clc;        % Clears the command window.
clear all;  % Clears all variables from the workspace.

% Input parameters for the printing process and geometry
H_of_layers = 0.3;         % Layer height in mm.
Square_Side_dimansion = 10; % Original side dimension of the square in mm.
Square_height = 10;        % Height of the square part in mm.
Elbow_Degree = 90;        % Angle of the elbow bend in degrees.
R_1 = 5;                 % Inner radius of the elbow in mm.
R_2 = 15;                % Outer radius of the elbow in mm.
offset = [1, 0, 0.6];   % Offset of the print bed origin in x, y, and z.
filament_diameter = 1.75; % Diameter of the filament in mm.
nozzle_diameter = 0.5;    % Diameter of the nozzle in mm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Step One: Optimize Layer Height for the Square ---
lower_bound = 0.7 * H_of_layers; % Lower bound for layer height optimization (70% of base layer height).
upper_bound = 1.3 * H_of_layers; % Upper bound for layer height optimization (130% of base layer height).

best_number = lower_bound; % Initialize the best layer height with the lower bound.
smallest_remainder = mod(Square_height, best_number); % Calculate the remainder of Square_height divided by the initial best number.

% Iterate through the range of possible layer heights to find the one that minimizes the remainder when dividing Square_height. This ensures a whole number of layers for the square.
for i = lower_bound:0.01:upper_bound
    remainder = mod(Square_height, i);
    if remainder < smallest_remainder
        best_number = i;
        smallest_remainder = remainder;
    end
end

H_of_layers_square = best_number; % Optimized layer height for the square.
Layers_Number_square = floor(Square_height / H_of_layers_square); % Number of layers for the square.
Square_Side = Square_Side_dimansion - (2 * H_of_layers_square); %Adjusted square side to account for layer height (important for proper dimensions).% Define the uni-coordinate system

% Position of the BED coordinate system within the uni-coordinate system

% Define the print bed origin
BED_origin = offset;
BED_origin_2 = BED_origin; % a copy the origin for later use.
 % distance of the first point from orgin of the coordinate
% Define the vertices of the cube part within the BED coordinate system

% pre-initial layer
  vertices =  [
                (BED_origin (1)),                 (BED_origin(2)),                  BED_origin(3);
                (BED_origin (1) - Square_Side),   (BED_origin (2)),                 BED_origin(3);
                (BED_origin (1) - Square_Side),   (BED_origin (2) + Square_Side),   BED_origin(3);
                (BED_origin (1)),                 (BED_origin (2) + Square_Side),   BED_origin(3)
                ];



% --- Step Two: Calculations for the Elbow Section ---
Circumference_min = (Elbow_Degree / 360) * 2 * R_1 * pi; % Minimum circumference of the elbow.
Circumference_max = (Elbow_Degree / 360) * 2 * R_2 * pi; % Maximum circumference of the elbow.

% Optimize layer height for the elbow
lower_bound = 0.7 * H_of_layers;
upper_bound = 1.3 * H_of_layers;

best_number = lower_bound;
smallest_remainder = mod(Circumference_min, best_number);

for i = lower_bound:0.01:upper_bound
    remainder = mod(Circumference_min, i);
    if remainder < smallest_remainder
        best_number = i;
        smallest_remainder = remainder;
    end
end

H_of_layers_Elbow = best_number; % Optimized layer height for the elbow.
Elbow_Size = Square_Side_dimansion - (2 * H_of_layers_Elbow); %Adjusted elbow size to account for layer height.
Gamma_Radian = H_of_layers_Elbow / (R_2); % Angle per layer in radians.
Gamma = Gamma_Radian * (180 / pi); % Angle per layer in degrees.


%%Extra layers;
% No layers; calculating the amounts of bug OR empty eara; it should be almost zero

% Calculate the number of layers and other related parameters for the elbow
Min_layers = Circumference_min / H_of_layers_Elbow;
Min_layers_number = floor(Min_layers);
Max_layers = Circumference_max / H_of_layers_Elbow;
Max_layers_number = floor(Max_layers);
mode = mod(Max_layers, Min_layers);

All_Layers = Max_layers_number;
final_layers_Number = 0;

if mode > H_of_layers_Elbow
    final_layers_Number = floor(mode / H_of_layers_Elbow);
end
% Circumference of sector of a layer = H_of_layers_Elbow
Circumference_final_layers = final_layers_Number * H_of_layers_Elbow;
% nolayers
no_layers = mod(Max_layers, H_of_layers_Elbow);
% degree of no layers area :
nolayers_degree = no_layers / (R_2);

All_Layers_Elbow = All_Layers + final_layers_Number;
Gamma_Total = All_Layers_Elbow * Gamma;

Layers_Number_2_Full = floor(Circumference_min/H_of_layers_Elbow); % full layers in elbow
% calculating final layers of the elbow_shell
Printed_Circumference_min = H_of_layers_Elbow * Layers_Number_2_Full;

Printed_Circumference_max = (Circumference_max * Printed_Circumference_min / Circumference_min);
%%% initial layers = final layers
final_Circumference =  (Circumference_max - Printed_Circumference_max);


% degree between layers

% rotation degree; it should be almost 90 degree or  Elbow_Degree
rotation_degree = (Elbow_Degree - nolayers_degree) / All_Layers_Elbow;
% incomplete_ layers numbering;
% number of incomplete spaces,
% spaces between full layers, the printing starts with incomplete layers
incomplete_spaces = Min_layers_number + round((final_layers_Number + 0.1)/(final_layers_Number + 0.3));
% number of incomplete_layers in each incomplete_spaces
incomplete_layers_number = (Max_layers_number - Min_layers_number)/incomplete_spaces;

    All_printing_layers = (2 * Layers_Number_square) + All_Layers_Elbow;
    AAA = cell(All_printing_layers, 4);

% Start of the printing simulation %%
%creating matrix that designate which layer must be rotate;
box (1: Layers_Number_square) = zeros;
box (Layers_Number_square + 1: Layers_Number_square+All_Layers_Elbow) = ones;
box (Layers_Number_square+All_Layers_Elbow+ 1: All_printing_layers ) = zeros;


j  = offset (3); % It must be always zero %% movement of Z in each step
Total_theta = 0; % It must be always zero
counting_1 = 0;
Save_matrix = [0;0;0];
small_theta_beteen = 0  ;  % small_theta_beteen rotated line and actual liNe
%% start the main loop
 for g  = 1:length(box); % number of layers %% length(box)
      AAA {g,1} = g;
      AAA {g,3} = 4; % number of line in incompelete layers
      AAA {g,4} = 0; %to find the rotating layers
  rotation_Degree = box (g)*rotation_degree; % rotation 0 or D-Of_rotation
  Total_theta = Total_theta + rotation_Degree;
%%% if situation for rotation;
% Define the translation to move x=1 axis to origin
translateToOrigin = [1 0 0 -offset(1); 0 1 0 0; 0 0 1 0; 0 0 0 1];

% Define the translation to move back after rotation
translateBack = [1 0 0 offset(1); 0 1 0 0; 0 0 1 0; 0 0 0 1];
if rotation_Degree > 0 % rotation bigger than 0
%  cla(axis1); %clear privious axis and its content
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

end
%%% End of if ;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the equation of a line  (a lyer which in going to be printed fully in the future)
% Calculate the slope (m)
% Ralpha = angle between curren layer to be printed with  the line  (a lyer which in going to be printed fully in the future)

% Elbow_printing
if box (g) == 1;
  counting_1 = counting_1 +1;


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

                x = 0;
                X = 0;
                Xt = Elbow_Size;
            else % Calculate x for other values of Ralpha
            x = (z -b)/m;
            X = x - R_1;
            Xt = Elbow_Size - X;
            AAA {g,3} = 3; %to not connect to end points
            end

BED_origin_2 (1) = BED_origin (1) + x_new;

% rotation and printing incomplete layers layers
vertices =  [
                (BED_origin_2 (1) - Xt),      (BED_origin(2)),                  j;
                (BED_origin_2 (1) ),          (BED_origin (2)),                 j;
                (BED_origin_2 (1) ),          (BED_origin (2) + Elbow_Size),   j;
                (BED_origin_2 (1) - Xt),      (BED_origin (2) + Elbow_Size),   j
                ];
%rotatedCubeVertices = (transformMatrix * [vertices, ones(size(vertices, 1), 1)]')';
 %          vertices = rotatedCubeVertices

         AAA {g,2} = vertices;

        % Close the square by connecting the last point to the first point


      pause(0.1)
         else
     j = j + H_of_layers_square;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% printing layers_square
if box (g) == 0;
   counting_22 = g;
   vertices =  [
              (BED_origin_2 (1)),                (BED_origin(2)),                  j;
              (BED_origin_2 (1) - Square_Side),  (BED_origin (2)),                 j;
              (BED_origin_2 (1) - Square_Side),  (BED_origin (2) + Square_Side),   j;
              (BED_origin_2 (1)),                (BED_origin (2) + Square_Side),   j
                ];
%rotatedCubeVertices = (transformMatrix * [vertices, ones(size(vertices, 1), 1)]')';
 %          vertices = rotatedCubeVertices

        % Close the square by connecting the last point to the first point
         AAA {g,2} = vertices;
end
end



%% start G-code writing:
%g_code_MATRIX =
for i = Layers_Number_square + 1:Layers_Number_square + All_Layers_Elbow;

 AAA{i, 4} = rotation_degree * (i - Layers_Number_square); % Rotation angle for each layer

end
fileID = fopen('generated_gcode.gcode', 'w');
% Write the initial setup commands
fprintf(fileID, ' ; Layer height=%.3f\n', H_of_layers);
fprintf(fileID, ' ; Modified Layer height for printing the first and last part, cubes = %.3f\n', H_of_layers_square);
fprintf(fileID, ' ; Modified Layer height for printing the second part, Elbow = %.3f\n', H_of_layers_Elbow);

fprintf(fileID, ' ; Modified the tube side size based on the Layer height = %.3f\n', Square_Side);
fprintf(fileID, ' ; Modified the elbow (R_2 - R_1) size based on the Layer height = %.3f\n', Elbow_Size);
fprintf(fileID, ' ; Angle of the Elbow = %.3f\n', Elbow_Degree);
fprintf(fileID, ' ; filament diameter = %.3f\n', filament_diameter);
fprintf(fileID, ' ; Angle of the Elbow = %.3f\n', Elbow_Degree);
fprintf(fileID, ' ; nozzle diameter = %.3f\n', nozzle_diameter);
fprintf(fileID, ' ; Modified Layer height for printing the first and last part, cubes = %.3f\n', H_of_layers_square );
fprintf(fileID, 'G21 ; Set units to millimeters\n');
fprintf(fileID, 'G90 ; Absolute positioning\n');
E = 0;
for layer = 1:All_printing_layers
    points = AAA{layer, 2}; % Get the points for the current layer
    rotation_angle = AAA{layer, 4}; % Get the rotation angle for the current layer
    fprintf(fileID, ';LAYER:%d\n', layer);
    % Move to the first point with rotation
    x = points(1, 1);
    y = points(1, 2);
    z = points(1, 3);
    fprintf(fileID, 'G0 F3600 X%.3f Y%.3f Z%.3f B%.2f\n', x, y, z, rotation_angle);

    % Loop through each point in the layer
    Point_index = 1: AAA{layer, 3};

% Perform the shift
    Point_index = [Point_index(2:end), Point_index(1)];
    for i = 1:length (Point_index)
    % number of G1 code; becuase the incomplete_layers are not complete square.
        point_idx =  Point_index (i);
        x = points(point_idx, 1);
        y = points(point_idx, 2);
        z = points(point_idx, 3);
        L = sqrt(x^2 + y^2); % Travel distance in XY plane
        if AAA{layer, 4} == 0
          layer_height = H_of_layers_square;
           else
           layer_height = H_of_layers_Elbow;
        end
           Ee = (4 * L * layer_height * nozzle_diameter) / (pi * (filament_diameter)^2);
            E = E + Ee;
        if point_idx == 2
            fprintf(fileID, 'G1 F600 X%.3f Y%.3f  E%.3f ; Move to point %d \n', x, y, E, point_idx);
        else
            fprintf(fileID, 'G1 X%.3f Y%.3f E%.3f ; Move to point %d \n', x, y, E, point_idx);
        end
    end



end

% Finish the G-code file
fprintf(fileID, 'M30 ; End of program\n');
fclose(fileID);

disp('G-code generation complete. Check the file generated_gcode.gcode.');




