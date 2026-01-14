#!/bin/bash



# Specify the folder path

folder_path="../rawdata"



# Specify the output file path

output_file="md5_checksums.txt"



# Change directory to the specified folder

cd "$folder_path" || exit



# Loop through each .tar.gz file

for file in *.fastq.gz; do

    # Calculate MD5 hash

    md5sum=$(md5sum "$file" | awk '{print $1}')



    # Append the file name and MD5 hash to the output file

    echo "$file $md5sum" >> "$output_file"

done



echo "MD5 checksums saved in $output_file"


