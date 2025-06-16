#!/usr/bin/env bash

###############################################################################
# Copyright (c) The JETSCAPE Collaboration, 2018
#
# For the list of contributors see AUTHORS.
#
# Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
#
# or via email to bugs.jetscape@gmail.com
#
# Distributed under the GNU General Public License 3.0 (GPLv3 or later).
# See COPYING for details.
##############################################################################

# clone the lbt-tables repository
repo_dir="LBT-tables-repo"
git clone https://github.com/JETSCAPE/LBT-tables.git "$repo_dir"

# temporary directory where the split LBT tables will be extracted
split_extract_dir="tmp/LBT-tables-split-extract"

# where the final unsplit LBT tables will be placed
unsplit_dir="LBT-tables"

if [ ! -d "$repo_dir/LBT-tables-split-archive" ]; then
    echo "Error archive directory does not exist."
    exit 1
fi

mkdir -p $split_extract_dir
mkdir -p $unsplit_dir

# extract all the .tar.gz files to the temporary directory
for file in "$repo_dir/LBT-tables-split-archive"/*.tar.gz; do
    if [ -f "$file" ]; then
        tar -xzf "$file" -C "$split_extract_dir"
        echo "Extracted '$file' to '$split_extract_dir'"
    else
        echo "Error '$file' does not conform to .tar.gz format."
        exit 1
    fi
done

# an associative array (map) for base names. eg. myFile.part0 and 
# myFile.part1 will share the base name myFile
declare -A base_names

# populate the base_names map with the unique base names
for file in "$split_extract_dir"/*; do
    filename=$(basename "$file")
    if [[ "$filename" =~ \.part([0-9]+)$ ]]; then
        base_name="${filename%.part*}"
        base_names["$base_name"]=1
    else
        echo "Error: '$filename' - all files must end with .partN, where N is a positive integer."
    fi
done

for base_name in "${!base_names[@]}"; do
    # sort parts numerically (not lexicographically) so a .part10 will not be appended before .part2
    part_files=($(ls "$split_extract_dir"/"$base_name".part* 2>/dev/null | sort -V))

    # if a file is small enough that it has only one part (.part0), copy it to the unsplit directory
    if [ ${#part_files[@]} -eq 1 ] && [[ "${part_files[0]}" == "$split_extract_dir/$base_name.part0" ]]; then
        cp "${part_files[0]}" "$unsplit_dir/$base_name"
        echo "Copied '${part_files[0]}' to '$unsplit_dir/$base_name'"
    else
        # concatenate the multiple parts, writing the file to the unsplit directory
        output_file="$unsplit_dir/$base_name"
        > "$output_file"
        for part in "${part_files[@]}"; do
            cat "$part" >> "$output_file"
            echo "Added '$part' to '$output_file'"
        done
    fi
done

# remove the temporary extracted split files
rm -rf tmp

# remove the local lbt-tables repository if not wanted
# rm -rf "$repo_dir"
