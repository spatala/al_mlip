#! /bin/bash

cnt=0

gb_folder="/Al_S17_0_N1_0_-1_4_N2_0_-1_-4"
full_gb_path="/home/leila/Desktop/symm_tilts/symm_tilts_100"
# gb_folder="Al_S5_0_N1_0_-1_3_N2_0_-1_-3"
# full_gb_path="/home/leila/Leila_sndhard/codes/GBMC-master"
for i in `ls -d $full_gb_path/$gb_folder/`;

do
   let cnt=cnt+1
   kk=$gb_folder"_"$(basename $i)
   echo $kk

   matlab -nodisplay -r "GBMC_RUN('$i'); quit"

   echo $i
done
