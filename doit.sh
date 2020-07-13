#!/bin/bash

input_dir="./Fakes/Swope/Galaxy_Fakes"

gal_bin_arr=(
"13.0_13.5"
"13.0_13.5"
"13.0_13.5"
"13.5_14.0"
"13.5_14.0"
"13.5_14.0"
"14.0_14.5"
"14.0_14.5"
"14.0_14.5"
"14.5_15.0"
"14.5_15.0"
"14.5_15.0"
"15.0_15.5"
"15.0_15.5"
"15.0_15.5"
"15.5_16.0"
"15.5_16.0"
"15.5_16.0"
"16.0_16.5"
"16.0_16.5"
"16.0_16.5"
"16.5_17.0"
"16.5_17.0"
"16.5_17.0"
"17.0_17.5"
"17.0_17.5"
"17.0_17.5"
)

bright_arr=(
18
20
21
18
20
21
18
20
21
18
20
21
18
20
21
18
20
21
18
20
21
18
20
21
18
20
21
)

dim_arr=(
21
22
25
21
22
25
21
22
25
21
22
25
21
22
25
21
22
25
21
22
25
21
22
25
21
22
25
)

fwhm_arr=(
3.0
3.0
3.0
3.0
3.0
3.0
3.0
3.0
3.0
3.0
3.0
3.0
3.0
3.0
3.0
3.0
3.0
3.0
3.0
3.0
3.0
3.0
3.0
3.0
3.0
3.0
3.0
)

global_start=$SECONDS
arr_len=${#gal_bin_arr[@]}
start_i=3 # use to skip to desired record

for i in "${!gal_bin_arr[@]}"; do # i==index, not object

  if [ $i -ge $start_i ]; then

    start=$SECONDS

    gal_bin=${gal_bin_arr[${i}]}
    gal_fake_bright=${bright_arr[${i}]}
    gal_fake_dim=${dim_arr[${i}]}
    fwhm_factor=${fwhm_arr[${i}]}

    iterations=$(expr $i + 1)
    iteration_start=${iterations}

    msg="Processing '${gal_bin}' between ${gal_fake_bright} and ${gal_fake_dim} with fwhm_multipier=${fwhm_factor}"
    msg="${msg} [${iterations}/${arr_len}] ..."
    echo "${msg}"

    # Do work...
    python ./DetEff_StaticFile.py \
    --stage plant,photpipe \
    --image_list ${input_dir}/${gal_bin}_images.txt \
    --template_list ${input_dir}/${gal_bin}_temps.txt \
    --iteration_start ${iteration_start} \
    --iterations ${iterations} \
    --gal_fake_mag_range ${gal_fake_bright} ${gal_fake_dim} 5000 0.2 \
    --gal_bin_to_process ${gal_bin} \
    --gal_fake_fwhm_factor ${fwhm_factor} \
    --plant_in_galaxies

    duration=$(( SECONDS - start ))
    msg="... done Processing '${gal_bin}' between ${gal_fake_bright} and ${gal_fake_dim}"
    msg="${msg} fwhm_multipier=${fwhm_factor}. [${iterations}/${arr_len}] elapsed: ${duration} sec."
    echo "${msg}"
    echo "" #newline

  else
    echo "Skipping index ${i} ..."
    echo "" #newline
  fi

done

global_duration=$(( SECONDS - start ))
echo "Full process elapsed: ${global_duration} sec."
