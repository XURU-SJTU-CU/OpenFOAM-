# OpenFOAM-Post Processing

## Some basics for postProcessing water level of OpenFOAM
the command line to excute the functions included in controlDict file
```
solverName  -postProcess > log.solverPost
```

## Post process functions

### IsoSurface postProcess
To compute the isosurface of a variable: for instance to compute the water surface isoSurface, one need to add the following to controlDict
```
functions
{

    surfaces1
    {
        type            surfaces;
        libs            (sampling);
        interpolationScheme cell;
        writeControl    writeTime;
        writeInterval   1;
        surfaceFormat   raw; //vtk, raw
        fields
        (
            alpha.water
        );
        surfaces
        {
            mySurface1
            {
                type            isoSurfaceCell;
                isoField        alpha.water;
                isoValues       (0.5);
                interpolate      true;
                regularise       true;
            }
        }
    }
}
```
then a surface data *.raw will be computed in postProcessing/surface1/0/  folder. If the file is exported to be vtk format, then it can be viewed in Paraview, while not directly accessible through visual-studio-code.

### Compute the water level of a specified position
The main idea to acquire the water level is that: find the closest position for example (x_0, y_0) in the *.raw file and output the z position accordingly, and minus the staring z position.
While if the case use dynamicMesh, then the x y positions will change with time if we want the water level of a certain point fixed in the moving frame. 

Say if the moving body moves according to a file which specifies the positions and rotations at each instant, 
The python code to obtain the water level specifying two points fixed in moving body frame Point1 and Point2 can be written as following:
```
import re
import os
import numpy as np
import matplotlib.pyplot as plt

def read_time_position(file_path, time_array, offset=(0, 0, 0)):
    """
    Reads a file containing time and position data, and retrieves interpolated positions for a specified time array.

    :param file_path: Path to the file containing time and position data.
    :param time_array: List or numpy array of target times to extract positions for.
    :param offset: Tuple of offsets (x_offset, y_offset, z_offset) to add to the output positions.
    :return: Dictionary with times as keys and corresponding interpolated positions (with offset) as values.
    """
    data = {}

    with open(file_path, 'r') as file:
        for line in file:
            # Match lines with time and position data using regex
            match = re.match(r"\(([^\s]+) \(\(([^\)]+)\) \(.*\)\)\)", line)
            if match:
                time = float(match.group(1))
                position = tuple(map(float, match.group(2).split()))

                # Store the time and position in the dictionary
                data[time] = position

    # Sort the data by time
    times = np.array(sorted(data.keys()))
    positions = np.array([data[t] for t in times])

    # Interpolate positions for the specified time array
    result = {}
    for target_time in time_array:
        if target_time < times[0] or target_time > times[-1]:
            result[target_time] = None  # Outside the interpolation range
        else:
            # Interpolate x, y, z separately
            interpolated_position = tuple(
                np.interp(target_time, times, positions[:, i]) + offset[i]
                for i in range(3)
            )
            result[target_time] = interpolated_position

    return result

def read_raw_file(file_path):
    """
    Reads a .raw file with isosurface data and returns the data as a numpy array.

    :param file_path: Path to the .raw file.
    :return: Numpy array of isosurface data.
    """
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            # Skip header lines
            if line.startswith('#') or line.strip() == "":
                continue
            # Parse the data lines
            values = list(map(float, line.split()))
            data.append(values)
    return np.array(data)

def find_nearest_z(data, x_target, y_target):
    """
    Finds the nearest z-value for a given (x, y) coordinate from the isosurface data.
    If there are multiple closest points, chooses the lowest z-value.

    :param data: Numpy array of isosurface data with columns [x, y, z, alpha.water].
    :param x_target: Target x-coordinate.
    :param y_target: Target y-coordinate.
    :return: The nearest z-value.
    """
    # Calculate the distance from the target (x, y) to all (x, y) in the data
    distances = np.sqrt((data[:, 0] - x_target)**2 + (data[:, 1] - y_target)**2)

    # Find the minimum distance
    min_distance = np.min(distances)

    # Get all indices with the minimum distance
    nearest_indices = np.where(distances == min_distance)[0]

    # Find the lowest z-value among the nearest points
    nearest_z = np.min(data[nearest_indices, 2])

    return nearest_z

def get_time_directories(base_path):
    """
    Reads all time directories from the base path and returns them as a sorted list of floats.

    :param base_path: Path containing time directories.
    :return: Sorted list of time directories as floats.
    """
    directories = []
    for entry in os.listdir(base_path):
        try:
            directories.append(float(entry))
        except ValueError:
            continue  # Skip non-numeric directories
    return sorted(directories)

# Combine both functionalities
file_path_time_position = "./constant/6DoF.dat"  # Replace with your file path
time_array = get_time_directories("./")  # Read times from existing directories

offset_point1 = (0, 0.1, 0)  # Example offset values for point1
offset_point2 = (0.1, 0.1, 0)  # Example offset values for point2

positions_point1 = read_time_position(file_path_time_position, time_array, offset_point1)
positions_point2 = read_time_position(file_path_time_position, time_array, offset_point2)

point1_zWaterLevel = []
point2_zWaterLevel = []

times = []

for time, position in positions_point1.items():
    if position:
        x, y, z = position
        directory_time = int(time) if time.is_integer() else time  # Use integer format if time is whole number
        file_path_iso = f"./postProcessing/surfaces1/{directory_time}/alpha.water_mySurface1.raw"  # Adjust file path based on time
        iso_data = read_raw_file(file_path_iso)
        nearest_z = find_nearest_z(iso_data, x, y)
        left_wall_z.append(nearest_z)
        times.append(time)

for time, position in positions_point2.items():
    if position:
        x, y, z = position
        directory_time = int(time) if time.is_integer() else time  # Use integer format if time is whole number
        file_path_iso = f"./postProcessing/surfaces1/{directory_time}/alpha.water_mySurface1.raw"  # Adjust file path based on time
        iso_data = read_raw_file(file_path_iso)
        nearest_z = find_nearest_z(iso_data, x, y)
        right_wall_z.append(nearest_z)

# Save data to a text file
output_file = "waterLevel_z_vs_time.txt"
with open(output_file, "w") as f:
    f.write("Time Point1 Point2 water level\n")
    for t, lw_z, rw_z in zip(times, point1_zWaterLevel, point2_zWaterLevel):
        f.write(f"{t:.6f} {lw_z:.6f} {rw_z:.6f}\n")

print(f"Data saved to {output_file}")

# Plotting the z-values
plt.figure(figsize=(10, 6))
plt.plot(times, point1_zWaterLevel, label="Point1 water level Z", marker="o")
plt.plot(times, point2_zWaterLevel, label="Point2 water level Z", marker="s")
plt.xlabel("Time[s]")
plt.ylabel("water level Z[m]")
plt.title("water level Z vs Time for Point1 and Point2")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("waterLevel_z_vs_time.png", dpi=300)
print("Plot saved as waterLevel_z_vs_time.png")
```
