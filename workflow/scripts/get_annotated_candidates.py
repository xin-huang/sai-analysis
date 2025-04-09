# Copyright 2025 Xin Huang
#
# GNU General Public License v3.0
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html


import pandas as pd
import os

# Input files from snakemake
outlier_file = snakemake.input.outliers  # Path to the outlier file
annotation_files = snakemake.input.annotation  # List of multianno files
output_file = snakemake.output.annotated_candidates  # Path to save results

with open(output_file, "w"):
    pass  # Ensure the output file is created

# Read the outlier file
outlier_df = pd.read_csv(outlier_file, sep="\t")

# Convert the annotation_files list into a dictionary for quick access
multianno_dict = {
    os.path.basename(f).split(".chr")[1].split(".")[0]: f for f in annotation_files
}

# Create a set to store unique (Chrom, Position) pairs
candidate_positions = set()

# Extract unique (Chrom, Position) pairs from U candidates and Q candidates
for _, row in outlier_df.iterrows():
    chrom = str(row["Chrom"]).removeprefix("chr")  # Convert chromosome to string for consistency

    # Extract candidate positions
    candidates = str(row["Overlapping Candidate"]).split(",")

    # Add (Chrom, Position) pairs to the set
    for c in candidates:
        pos = int(c.split(":")[1])
        candidate_positions.add((chrom, pos))

# List to store filtered variants
filtered_variants = []

# Process each unique chromosome
for chrom in set(chrom for chrom, _ in candidate_positions):
    print(f"Processing chromosome {chrom}...")

    # Check if the chromosome has a corresponding annotation file
    if chrom not in multianno_dict:
        print(f"No annotation file for chromosome: {chrom}")
        continue

    # Load the multianno file
    multianno_file = multianno_dict[chrom]
    print(f"Reading annotation file: {multianno_file}")

    multianno_df = pd.read_csv(multianno_file, sep="\t")

    # Ensure the relevant columns are in the correct type
    multianno_df["Chr"] = multianno_df["Chr"].astype(str).apply(lambda x: x.removeprefix("chr"))
    multianno_df["Start"] = multianno_df["Start"].astype(int)
    multianno_df["End"] = multianno_df["End"].astype(int)

    # Get all positions relevant to this chromosome
    positions_for_chrom = {pos for c, pos in candidate_positions if c == chrom}

    # Filter variants where 'Start' matches any candidate position
    variants_in_region = multianno_df[multianno_df["Start"].isin(positions_for_chrom)]

    print(f"Filtered {len(variants_in_region)} variants for chromosome {chrom}.")

    # Append to the results list
    if not variants_in_region.empty:
        filtered_variants.append(variants_in_region)

# Combine all filtered variants into a single DataFrame
if filtered_variants:
    result_df = pd.concat(filtered_variants, ignore_index=True)
    result_df["Chr"] = result_df["Chr"].astype(int)
    result_df = result_df.sort_values(
        by=["Chr", "Start", "End"], ascending=[True, True, True]
    )
    # Save the filtered variants to a file
    result_df.to_csv(output_file, sep="\t", index=False)
