#!/bin/bash

input_dir="./Fakes/Swope/Galaxy_Fakes"
#gal_bin="13.0_13.5"
#gal_fake_bright=18
#gal_fake_dim=21
#fwhm_factor=3.0

#echo "Running gal fakes for (13.0, 13.5), gal_fake_mag_range=(${gal_fake_bright}, ${gal_fake_dim}) ..."
#python ./DetEff_StaticFile.py \
#  --stage plant,photpipe \
#  --image_list ${input_dir}/${gal_bin}_images.txt \
#  --template_list ${input_dir}/${gal_bin}_temps.txt \
#  --iteration_start 1 \
#  --iterations 1 \
#  --gal_fake_mag_range ${gal_fake_bright} ${gal_fake_dim} 5000 0.2 \
#  --gal_bin_to_process ${gal_bin} \
#  --gal_fake_fwhm_factor ${fwhm_factor} \
#  --plant_in_galaxies
#echo "... Done running gal fakes for (13.0, 13.5), gal_fake_mag_range=({$gal_fake_bright}, {$gal_fake_dim})"

# fake_params for each run:
# gal_bin, gal_fake_bright, gal_fake_dim, fwhm_factor
#fake_params=(
#("13.0_13.5", 18, 21, 3.0)
#("13.0_13.5", 20, 22, 3.0)
#)

gal_bin_arr=(
"13.0_13.5"
"13.0_13.5"
)

bright_arr=(
18
20
)

dim_arr=(
21
22
)

fwhm_arr=(
3.0
3.0
)

for i in `seq 0 ${#gal_bin_arr[@]}`; do

  start=$SECONDS

  gal_bin=${gal_bin_arr[${i}]}
  gal_fake_bright=${bright_arr[${i}]}
  gal_fake_dim=${dim_arr[${i}]}
  fwhm_factor=${fwhm_arr[${i}]}
  msg="Processing '${gal_bin}' between ${gal_fake_bright} and ${gal_fake_dim} with fwhm_multipier=${fwhm_factor} ..."
  echo "${msg}"

  # Do work...

  duration=$(( SECONDS - start ))
  msg="... done Processing '${gal_bin}' between ${gal_fake_bright} and ${gal_fake_dim}"
  msg="${msg} fwhm_multipier=${fwhm_factor}. Elapsed: ${duration} sec."
  echo "${msg}"
  echo "" #newline

done

#gal_fake_bright=20
#gal_fake_dim=22
#python ./DetEff_StaticFile.py \
#  --stage plant,photpipe \
#  --image_list ${input_dir}/${gal_bin}_images.txt \
#  --template_list ${input_dir}/${gal_bin}_temps.txt \
#  --iteration_start 2 \
#  --iterations 2 \
#  --gal_fake_mag_range ${gal_fake_bright} ${gal_fake_dim} 5000 0.2 \
#  --gal_bin_to_process ${gal_bin} \
#  --gal_fake_fwhm_factor ${fwhm_factor} \
#  --plant_in_galaxies
#echo "... Done running gal fakes for (13.0, 13.5), gal_fake_mag_range=({$gal_fake_bright}, {$gal_fake_dim})"
#
#gal_fake_bright=21
#gal_fake_dim=25
#python ./DetEff_StaticFile.py \
#  --stage plant,photpipe \
#  --image_list ${input_dir}/${gal_bin}_images.txt \
#  --template_list ${input_dir}/${gal_bin}_temps.txt \
#  --iteration_start 3 \
#  --iterations 3 \
#  --gal_fake_mag_range ${gal_fake_bright} ${gal_fake_dim} 5000 0.2 \
#  --gal_bin_to_process ${gal_bin} \
#  --gal_fake_fwhm_factor ${fwhm_factor} \
#  --plant_in_galaxies
#echo "... Done running gal fakes for (13.0, 13.5), gal_fake_mag_range=({$gal_fake_bright}, {$gal_fake_dim})"
