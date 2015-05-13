import sys

# Open both fastq files and make single array with combined reads
with open(sys.argv[1]) as f1:
    file1 = f1.readlines()

with open(sys.argv[2]) as f2:
    file2 = f2.readlines()

combined_reads = [a.rstrip()+" "+b for (a,b) in zip(file1, file2)]

# Go through reads and record which reads to remove
seq_count = 1
line_count = 1
remove = []
for line in combined_reads:
    if seq_count == 5:
		seq_count = 1
    if seq_count == 2 and sys.argv[3] in line:
		remove.append(line_count - 1)
		remove.append(line_count)
		remove.append(line_count + 1)
		remove.append(line_count + 2)
    seq_count += 1
    line_count += 1


# Write new files, removing unwanted reads
chopped_file1 = open(sys.argv[4], "w")
line_count = 1
for line in file1:
    if line_count not in remove:
        chopped_file1.write(line)
    line_count = line_count + 1

chopped_file2 = open(sys.argv[5], "w")
line_count = 1
for line in file2:
    if line_count not in remove:
        chopped_file2.write(line)
    line_count = line_count + 1

