import os, sys
import matplotlib.pyplot as plt

# Function to count lines containing ">" in a file
def count_lines_with_char(file_path, char=">"):
    count = 0
    with open(file_path, 'r') as file:
        for line in file:
            if char in line:
                count += 1
    return count

# Directory containing the files (adjust the path as needed)
directory = sys.argv[1]

# List to store file numbers and corresponding counts
file_numbers = []
line_counts = []

# Iterate over each file in the directory
for filename in os.listdir(directory):
    if filename.endswith(".fasta"):  # Adjust the extension as needed
        # Extract the number from the filename (assuming filenames are in format like file1.txt)
        file_number = int(''.join(filter(str.isdigit, filename)))
        
        # Full path to the file
        file_path = os.path.join(directory, filename)
        
        # Count lines containing ">"
        count = count_lines_with_char(file_path)
        
        # Append the results to the lists
        file_numbers.append(file_number)
        line_counts.append(count)

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(file_numbers, line_counts, marker='o')
plt.xlabel('File Number')
plt.ylabel('Number of Lines Containing ">"')
plt.title('Number of Lines Containing ">" in Each File')
plt.grid(True)

# Save the plot to a file
output_file = 'plot_90_100_1.png'  # Change the file name and extension as needed
plt.savefig(output_file)

print(f"Plot has been saved to '{output_file}'")
