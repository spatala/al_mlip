#! /bin/bash

cnt=0

full_gb_path="/home/leila/Desktop/symm_tilts/symm_tilts_100"
regex="Al_S([0-9]+)_"

for i in `ls -d $full_gb_path/*`;
do
	if [[ $i =~ $regex ]]
	then
		sig_num="${BASH_REMATCH[1]}"
		if   [ "$sig_num" -le 200 ]
		then
			echo $sig_num
			let cnt=cnt+1
		    kk=$gb_folder"_"$(basename $i)
		    echo $kk

		    matlab -nodisplay -r "GBMC_RUN('$i'); quit"

		    echo $i
		fi
	fi
done
