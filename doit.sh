#!/bin/bash


# Swope
#input_dir="./Fakes/Swope/Galaxy_Fakes"
#root_path="/data/LCO/Swope/workstch"
#field_name_start="s005"
#log_base="logstch"
#work_base="workstch"

# Thacher
input_dir="./Fakes/Thacher/Galaxy_Fakes"
root_path="/data2/THACHER/workspace"
field_name_start="t005"
log_base="logs"
work_base="workspace"

#gal_bin_arr=(
#"13.0_13.5" # 0
#"13.0_13.5" # 1
#"13.0_13.5" # 2
#"13.5_14.0" # 3
#"13.5_14.0" # 4
#"13.5_14.0" # 5
#"14.0_14.5" # 6
#"14.0_14.5" # 7
#"14.0_14.5" # 8
#"14.5_15.0" # 9
#"14.5_15.0" # 10
#"14.5_15.0" # 11
#"15.0_15.5" # 12
#"15.0_15.5" # 13
#"15.0_15.5" # 14
#"15.5_16.0" # 15
#"15.5_16.0" # 16
#"15.5_16.0" # 17
#"16.0_16.5" # 18
#"16.0_16.5" # 19
#"16.0_16.5" # 20
#"16.5_17.0" # 21
#"16.5_17.0" # 22
#"16.5_17.0" # 23
#"17.0_17.5" # 24
#"17.0_17.5" # 25
#"17.0_17.5" # 26
#)

#bright_arr=(
#18 # 0
#20 # 1
#21 # 2
#18 # 3
#20 # 4
#21 # 5
#18 # 6
#20 # 7
#21 # 8
#18 # 9
#20 # 10
#21 # 11
#18 # 12
#20 # 13
#21 # 14
#18 # 15
#20 # 16
#21 # 17
#18 # 18
#20 # 19
#21 # 20
#18 # 21
#20 # 22
#21 # 23
#18 # 24
#20 # 25
#21 # 26
#)
#
#dim_arr=(
#21 # 0
#22 # 1
#25 # 2
#21 # 3
#22 # 4
#25 # 5
#21 # 6
#22 # 7
#25 # 8
#21 # 9
#22 # 10
#25 # 11
#21 # 12
#22 # 13
#25 # 14
#21 # 15
#22 # 16
#25 # 17
#21 # 18
#22 # 19
#25 # 20
#21 # 21
#22 # 22
#25 # 23
#21 # 24
#22 # 25
#25 # 26
#)

#fwhm_arr=(
#3.0 # 0
#3.0 # 1
#3.0 # 2
#3.0 # 3
#3.0 # 4
#3.0 # 5
#3.0 # 6
#3.0 # 7
#3.0 # 8
#3.0 # 9
#3.0 # 10
#3.0 # 11
#5.0 # 12
#5.0 # 13
#5.0 # 14
#3.0 # 15
#3.0 # 16
#3.0 # 17
#3.0 # 18
#3.0 # 19
#3.0 # 20
#3.0 # 21
#3.0 # 22
#3.0 # 23
#3.0 # 24
#3.0 # 25
#3.0 # 26
#)

global_start=$SECONDS
#arr_len=${#gal_bin_arr[@]}
#start_i=16 # index to start (i.e. skips indices < start_i)
#stop_i=16 # index to stop after (i.e. skips indices > stop_i)

# Hack for single run use...
#dir_append=4

START=2
END=10

#for i in $( seq $START $END ); do # i==index, not object
#for i in "${!gal_bin_arr[@]}"; do # i==index, not object

#  if [ $i -lt $start_i -o $i -gt $stop_i ]; then
#    echo "Skipping index ${i} ..."
#    echo "" #newline
#    continue
#  else

start=$SECONDS

#  gal_bin=${gal_bin_arr[${i}]}
#  gal_fake_bright=${bright_arr[${i}]}
#  gal_fake_dim=${dim_arr[${i}]}
#  fwhm_factor=${fwhm_arr[${i}]}
gal_bin="13.0_13.5"
gal_fake_bright="18"
#gal_fake_dim="23"
gal_fake_dim="22"
fwhm_factor="5"

convolve_which='b'
#    iterations=$(($i + 1))
#    iteration_start=${iterations}
#  iterations=${dir_append}
#  iteration_start=${dir_append}
#  current_dir=$i
#  end_dir=$i



msg="Processing '${gal_bin}' between ${gal_fake_bright} and ${gal_fake_dim} with fwhm_multipier=${fwhm_factor}"
#  msg="${msg} [${iterations}/${arr_len}] ..."
msg="${msg} [${iterations}/${END}] ..."
echo "${msg}"
#
# Do work... #  --stage plant,photpipe \
python ./DetEff_StaticFile.py \
--stage plant,photpipe \
--root_path ${root_path} \
--log_base ${log_base} \
--work_base ${work_base} \
--field_name_start ${field_name_start} \
--image_list ${input_dir}/${gal_bin}_images.txt \
--template_list ${input_dir}/${gal_bin}_temps.txt \
--subdir_start ${START} \
--subdir_end ${END} \
--convolve_which ${convolve_which} \
--gal_fake_mag_range ${gal_fake_bright} ${gal_fake_dim} 5000 0.2 \
--gal_bin_to_process ${gal_bin} \
--fake_fwhm_factor ${fwhm_factor} \
--plant_in_galaxies

duration=$(( SECONDS - start ))
msg="... done Processing '${gal_bin}' between ${gal_fake_bright} and ${gal_fake_dim}"
#  msg="${msg} fwhm_multipier=${fwhm_factor}. [${iterations - }/${arr_len}] elapsed: ${duration} sec."
echo "${msg}"
echo "" #newline
#
#  ((dir_append++))

#  fi

#done

global_duration=$(( SECONDS - start ))
echo "Full process elapsed: ${global_duration} sec."
