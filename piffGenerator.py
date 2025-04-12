import numpy as np
import os

def generate_pif(input_path, file_name):
    
    # File paths
    
    input_file_path = os.path.join(input_path, file_name)
    
    oinput_parent = os.path.dirname(input_path)

    output_pif_path =os.path.join(oinput_parent, 'output.piff')
    
    # Load the input file and convert it into a 2D NumPy array
    with open(input_file_path, 'r') as file:
        data_matrix = np.array([list(map(int, line.split())) for line in file])

    # Get the dimensions of the data matrix
    xsize, ysize = data_matrix.shape

    def find_neighbors(matrix, x, y, visited):
        """ Helper function to find all connected cells of the same type """
        cell_value = matrix[x, y]
        stack = [(x, y)]
        cell_coordinates = []

        while stack:
            cx, cy = stack.pop()
            if (cx, cy) not in visited and matrix[cx, cy] == cell_value:
                visited.add((cx, cy))
                cell_coordinates.append((cx, cy))
                # Check the neighbors (up, down, left, right)
                if cx > 0:
                    stack.append((cx - 1, cy))
                if cx < xsize - 1:
                    stack.append((cx + 1, cy))
                if cy > 0:
                    stack.append((cx, cy - 1))
                if cy < ysize - 1:
                    stack.append((cx, cy + 1))

        return cell_coordinates

    # Track visited cells to avoid redundancy
    visited = set()
    cell_id = 1

    # Open the PIF file for writing
    with open(output_pif_path, 'w') as pif_file:
        pif_file.write(f"0 Medium 0 {xsize-1} 0 {ysize-1} 0 0\n")
        # Iterate through the data matrix to find unique cells
        for x in range(xsize):
            for y in range(ysize):
                if (x, y) not in visited and data_matrix[x, y] != 1:  # Skip the Medium cells
                    cell_coordinates = find_neighbors(data_matrix, x, y, visited)
                    min_x = min(cell_coordinates, key=lambda t: t[0])[0]
                    max_x = max(cell_coordinates, key=lambda t: t[0])[0]
                    min_y = min(cell_coordinates, key=lambda t: t[1])[1]
                    max_y = max(cell_coordinates, key=lambda t: t[1])[1]
                    # Write the cell information in PIF format
                    pif_file.write(f"{cell_id} CELL {min_x} {max_x} {min_y} {max_y} 0 0\n")
                    cell_id += 1

    # Confirmation message
    print(f"PIF file created successfully at {output_pif_path}")
