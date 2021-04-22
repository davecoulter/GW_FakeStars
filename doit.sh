#!/bin/bash

start=$SECONDS

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

START=1
END=10

gal_bin="13.5_14.0"
gal_fake_bright="18"
gal_fake_dim="22"
fwhm_factor="3"
convolve_which='b'

msg="Processing '${gal_bin}' between ${gal_fake_bright} and ${gal_fake_dim} with fwhm_multipier=${fwhm_factor}"
echo "${msg}"

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
echo "${msg}"
echo "" #newline

