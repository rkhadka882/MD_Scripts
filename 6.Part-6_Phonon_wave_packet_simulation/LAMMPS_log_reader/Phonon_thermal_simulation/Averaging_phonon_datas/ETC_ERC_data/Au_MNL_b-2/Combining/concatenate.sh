#!/bin/bash

# Read the job.dat file
while IFS=$'\t' read -r jid1 jid2
do
  # Combine files and rename
  file1="Energy_etc_erc_$jid1.dat"
  file2="Energy_etc_erc_$jid2.dat"
  combined="combined.txt"

  tail -n +9 "$file2" | cat "$file1" - > "$combined"
  mv "$combined" "$file1"

done < "job.dat"

