# G-code-generator and Printing simulator of a 5-axis printer
![octave-gui_kppVYGOOGi](https://github.com/user-attachments/assets/650b04e2-e81e-41a7-a59f-f4d4f625abdb)

## 5-Axis Slicing for FDM 3D Printers (GNU Octave/MATLAB) - G-Code Generation
![Screenshot201632](https://github.com/user-attachments/assets/95ffe556-f8ff-4060-8a9c-9c90e08f55b3)

This project implements a simplified 5-axis slicing algorithm and generates G-code for a Cartesian FDM 3D printer with simulated B-axis (rotation around the Y-axis) capabilities. It's designed to allow printing of curved beams or similar geometries without the need for support structures. There are several solution to do so This code focus on solving the program based on constant layers' height method. 

## Project Overview

This code addresses the challenge of printing curved parts on a standard 3-axis FDM printer by simulating a 5-axis printing process. By rotating the printing plane during the print, the need for support structures can be eliminated for certain geometries. This project generates G-code that includes these rotational movements.
![image](https://github.com/user-attachments/assets/2fb1fb38-16f0-4dbf-99fc-1309ef098749)

Reference of the Picture: https://doi.org/10.48550/arXiv.2202.11426
## Features

*   **Simulated 5-Axis Slicing:** Simulates the slicing of a part and generates toolpaths with rotational adjustments around the Y-axis (B-axis).
*   **Layer Height Optimization:** Optimizes layer height for both the straight and curved sections to ensure a consistent number of layers and minimize rounding errors.
*   **G-Code Generation:** Generates G-code that includes X, Y, Z movements and B-axis rotations.
*   **Elbow Geometry Handling:** Specifically designed to handle a geometry consisting of straight sections connected by a 90-degree elbow.
*   **Adjusted Dimensions:** Adjusts the dimensions of the square and elbow sections based on the layer height to ensure accurate final dimensions.

## How It Works

The code takes the following inputs:

*   `H_of_layers`: Base layer height.
*   `Square_Side_dimansion`: Original side dimension of the square.
*   `Square_height`: Height of the square sections.
*   `Elbow_Degree`: Angle of the elbow (90 degrees in this example).
*   `R_1`: Inner radius of the elbow.
*   `R_2`: Outer radius of the elbow.
*   `offset`: Initial offset of the print bed.
*   `filament_diameter`: Diameter of the filament.
*   `nozzle_diameter`: Diameter of the nozzle.

The process is as follows:

1.  **Layer Height Optimization:** Optimizes the layer height for the square and elbow sections.
2.  **Geometry Calculations:** Calculates geometric parameters for the elbow, such as circumferences, angles, and the number of layers. The dimensions of the square and elbow are adjusted to account for the layer height.
3.  **Layer Generation and Rotation:** Generates the vertices for each layer and calculates the rotation angle for the elbow layers.
4.  **G-Code Generation:** Creates a `.gcode` file containing the toolpaths, including linear movements (G1), rapid moves (G0), and B-axis rotations. The code also calculates and includes extrusion amounts (E values).

## Getting Started
# To simulate the process of printing: 
1.  **Prerequisites:** GNU Octave (or MATLAB).
2.  **Clone the Repository:**
3. **Run the Code:** Open `printing_process_visualizer.m` in MATLAB or Octave and run it.


4.  **Run the Code:** Open the `G-code-generator.m` file in Octave/MATLAB and run it. This will generate a file named `generated_gcode.gcode` in the same directory.

## Code Structure

*   The code is contained within two `.m` files.
*   The code is structured into logical sections: input parameters, layer height optimization, elbow calculations, layer generation, and G-code writing.

## G-Code Output

The generated G-code includes:

*   **Header:** Contains comments with input parameters and setup commands (G21, G90).
*   **Layer Movements:** G0 and G1 commands for positioning the nozzle in X, Y, and Z.
*   **B-Axis Rotations:** B commands for rotating the printing plane.
*   **Extrusion:** E values for controlling the extruder.
*   **Layer Comments:** Comments indicating the start of each layer.


## Important Considerations

*   **Elbow Approximation:** The elbow is approximated by straight lines.
*   **Machine Compatibility:** The generated G-code may require adjustments depending on the specific 5-axis machine configuration. This code assumes a B-axis.
*   **No Collision Detection:** There is no collision detection implemented in this code.

![Screenshot192840](https://github.com/user-attachments/assets/bcc9583d-19ba-4ca4-bdbd-cfa8d064f288)

## Contributing
## Output results: 
![Screenshot 2024-12-10 233459](https://github.com/user-attachments/assets/02d91433-f7d8-4e54-a7c8-9554d0742ad0)
![Screenshot 2024-12-10 232424](https://github.com/user-attachments/assets/fe1de620-b66c-49ba-b585-648c2d748086)
Contributions are welcome! Please open an issue or submit a pull request.


