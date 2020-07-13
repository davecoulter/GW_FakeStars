#!/bin/bash

input_dir="./Fakes/Swope/Galaxy_Fakes"
gal_bin="13.0_13.5"
gal_fake_bright="18"
gal_fake_dim="21"

echo "Running gal fakes for (13.0, 13.5), gal_fake_mag_range=($gal_fake_bright, $gal_fake_dim) ..."
python ./DetEff_StaticFile.py \
  --stage plant,photpipe \
  --image_list $input_dir/$gal_bin_images.txt \
  --template_list $input_dir/$gal_bin_temps.txt \
  --iteration_start 1 \
  --iterations 1 \
  --gal_fake_mag_range $gal_fake_bright $gal_fake_dim 5000 0.2 \
  --gal_bin_to_process $gal_bin \
  --gal_fake_fwhm_factor 3.0 \
  --plant_in_galaxies
echo "... Done running gal fakes for (13.0, 13.5), gal_fake_mag_range=($gal_fake_bright, $gal_fake_dim)"

gal_fake_bright="20"
gal_fake_dim="22"
echo "Running gal fakes for (13.0, 13.5), gal_fake_mag_range=($gal_fake_bright, $gal_fake_dim) ..."
python ./DetEff_StaticFile.py \
  --stage plant,photpipe \
  --image_list $input_dir/$gal_bin_images.txt \
  --template_list $input_dir/$gal_bin_temps.txt \
  --iteration_start 1 \
  --iterations 1 \
  --gal_fake_mag_range $gal_fake_bright $gal_fake_dim 5000 0.2 \
  --gal_bin_to_process $gal_bin \
  --gal_fake_fwhm_factor 3.0 \
  --plant_in_galaxies
echo "... Done running gal fakes for (13.0, 13.5), gal_fake_mag_range=($gal_fake_bright, $gal_fake_dim)"

gal_fake_bright="21"
gal_fake_dim="25"
echo "Running gal fakes for (13.0, 13.5), gal_fake_mag_range=($gal_fake_bright, $gal_fake_dim) ..."
python ./DetEff_StaticFile.py \
  --stage plant,photpipe \
  --image_list $input_dir/$gal_bin_images.txt \
  --template_list $input_dir/$gal_bin_temps.txt \
  --iteration_start 1 \
  --iterations 1 \
  --gal_fake_mag_range $gal_fake_bright $gal_fake_dim 5000 0.2 \
  --gal_bin_to_process $gal_bin \
  --gal_fake_fwhm_factor 3.0 \
  --plant_in_galaxies
echo "... Done running gal fakes for (13.0, 13.5), gal_fake_mag_range=($gal_fake_bright, $gal_fake_dim)"
