import os

# Path to the directory containing the subdirectories file1, file2, ..., file40
base_dir = os.getcwd()

# Name of the file to read from each subdirectory
filename = 'BERRYCURV_KUBO.dat'

# Output file to store the combined data
output_file = 'BERRYCURV_KUBO-tot.dat'

# Initialize a list to store combined data
combined_data = []

# Loop through directories file1, file2, ..., file40
for i in range(1, 41):
    dir_name = f'file{i}'
    file_path = os.path.join(base_dir, dir_name, filename)
    
    # Check if the file exists
    if os.path.isfile(file_path):
        with open(file_path, 'r') as file:
            lines = file.readlines()
            # Append lines starting from the 26th line (index 25)
            combined_data.extend(lines[25:])

# Write combined data to the output file
with open(output_file, 'w') as outfile:
    outfile.writelines(combined_data)

print(f"Data combined into {output_file}")

