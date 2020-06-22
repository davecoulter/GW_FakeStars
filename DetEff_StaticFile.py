#!/usr/bin/env python
# D. Jones - 5/3/18
# D. Coulter - 6/1/18
from __future__ import print_function
import shutil
import time
import random

"""Get swope/DOPHOT detection efficiencies as a 
function of mag using PHOTPIPE"""

import glob
import numpy as np
from astropy.io import fits
import os
from txtobj import txtobj
import re
from runsex import runsex
import sys
from astLib import astCoords
import astropy.table as at

class FileAssociation():
    def __init__(self, image_file,
                 image_dcmp_file,
                 image_mask_file,
                 image_noise_file,
                 fake_image_file,
                 fake_image_dcmp_file,
                 fake_image_mask_file,
                 fake_image_noise_file,
                 template_file,
                 template_dcmp_file,
                 diff_dcmp_file,
                 segmap_file,
                 sexcat_good = None,
                 pixel_good = None):

        self.image_file = image_file
        self.image_dcmp_file = image_dcmp_file
        self.image_mask_file = image_mask_file
        self.image_noise_file = image_noise_file
        self.fake_image_file = fake_image_file
        self.fake_image_dcmp_file = fake_image_dcmp_file
        self.fake_image_mask_file = fake_image_mask_file
        self.fake_image_noise_file = fake_image_noise_file
        self.template_file = template_file
        self.template_dcmp_file = template_dcmp_file
        self.diff_dcmp_file = diff_dcmp_file
        self.segmap_file = segmap_file

        self.sexcat_good = sexcat_good
        self.pixel_good = pixel_good

    def __repr__(self):
        return str(self.__dict__)

class DetermineEfficiencies():

    def __init__(self, root_path, image_dir, template_dir, image_list, template_list):

        self.root_path = root_path
        self.image_dir = image_dir
        self.template_dir = template_dir

        self.image_path = "{0}/{1}/1".format(self.root_path, self.image_dir)
        self.template_path = "{0}/{1}/1".format(self.root_path, self.template_dir)
        self.log_path = self.image_path.replace("workstch", "logstch")

        self.image_list = image_list
        self.template_list = template_list

        # Unpack image list
        self.image_names = []
        with open(image_list, 'r') as input_file:
            self.image_names += input_file.read().splitlines()

        self.image_files = []
        for i in self.image_names:
            self.image_files.append("{0}/{1}".format(self.image_path, i))

        # Unpack image list
        self.template_names = []
        with open(template_list, 'r') as input_file:
            self.template_names += input_file.read().splitlines()

        self.template_files = []
        for t in self.template_names:
            self.template_files.append("{0}/{1}".format(self.template_path, t))

    def initialize(self, iteration, gal_bin_to_process):

        if self.options.plant_in_galaxies:
            galstr = 'gal_' + gal_bin_to_process

            tokens = gal_bin_to_process.split("_")
            self.start_mag = float(tokens[0])
            self.end_mag = float(tokens[1])

        else:
            galstr = ''

        # These directories depend on the iteration #
        if len(galstr) > 0:
            self.fake_image_dir = "{0}_fake_{1}_{2}".format(self.image_dir, galstr, iteration)
        else:
            self.fake_image_dir = "{0}_fake_{1}{2}".format(self.image_dir, galstr, iteration)

        self.fake_image_path = "{0}/{1}/1".format(self.root_path, self.fake_image_dir)
        self.fake_log_path = self.fake_image_path.replace("workstch", "logstch")
        self.fake_stitched_path = self.fake_log_path.replace("logstch", "stitched")

        if len(galstr) > 0:
            self.diff_dir_name = "%s_fake_%s_%s_%s" % (self.image_dir, galstr, iteration, self.template_dir)
        else:
            self.diff_dir_name = "%s_fake_%s%s_%s" % (self.image_dir, galstr, iteration, self.template_dir)

        self.diff_dir_path = "%s/%s/1" % (self.root_path, self.diff_dir_name)

        self.fake_mag_file = self.fake_log_path + "/fakemags.txt"
        self.out_match_file = self.fake_log_path + "/fakemeasuredmags.txt"
        self.efficiency_path = self.fake_log_path + "/outeff/"
        self.full_efficiency_path = self.diff_dir_path + "/outeff/"

        self.initialize_dirs_and_logs()
        self.file_associations = self.build_file_associations()

    def initialize_dirs_and_logs(self):

        if not os.path.exists(self.fake_image_path):
            os.makedirs(self.fake_image_path)

        if not os.path.exists(self.fake_log_path):
            os.makedirs(self.fake_log_path)

        if not os.path.exists(self.full_efficiency_path):
            os.makedirs(self.full_efficiency_path)

        if not os.path.exists(self.fake_stitched_path):
            os.makedirs(self.fake_stitched_path)

        # Create file to hold fake magnitudes
        if not os.path.isfile(self.fake_mag_file) or self.options.stage == 'plant':
            with open(self.fake_mag_file, 'w') as fout:
                print('# dcmpfile x y mag mag_gal mag_gal_GLADE GLADE_gal_id ', file=fout)

        if not os.path.isfile(self.out_match_file) or self.options.stage == 'getphot':
            with open(self.out_match_file, 'w') as fout:
                print('# dcmpfile x y mag_sim mag_gal detmag detmagerr', file=fout)

        if not os.path.exists(self.efficiency_path):
            os.makedirs(self.efficiency_path)

    def build_file_associations(self):

        file_associations = {}

        image_fields = [i.split('.')[0] for i in self.image_names]
        template_fields = [t.split('.')[0] for t in self.template_names]

        for i, image_name in enumerate(self.image_names):

            image_file = self.image_files[i]
            image_dcmp_file = image_file.replace(".fits", ".dcmp")
            if os.path.exists(image_file.replace('.fits', '.mask.fits.gz')):
                image_mask_file = image_file.replace('.fits', '.mask.fits.gz')
            else:
                image_mask_file = image_file.replace('.fits', '.mask.fits')
            if os.path.exists(image_file.replace('.fits', '.noise.fits.gz')):
                image_noise_file = image_file.replace('.fits', '.noise.fits.gz')
            else:
                image_noise_file = image_file.replace('.fits', '.noise.fits')
            try:
                template_index = np.where(np.asarray(template_fields) == image_fields[i])[0][0]
            except:
                continue
            template_name = self.template_names[template_index]
            template_file = self.template_files[template_index]
            template_dcmp_file = template_file.replace("fits", "dcmp")

            date_match = re.findall(r"ut1\d{5}", image_name)[0]
            fake_image_name = image_name.replace(date_match, "{0}_fake".format(date_match))
            fake_image_file = "%s/%s" % (self.fake_image_path, fake_image_name)
            fake_image_dcmp_file = "%s/%s" % (self.fake_image_path, fake_image_name.replace(".fits", ".dcmp"))
            fake_image_mask_file = "%s/%s" % (self.fake_image_path, fake_image_name.replace('.fits', '.mask.fits.gz'))
            fake_image_noise_file = "%s/%s" % (self.fake_image_path, fake_image_name.replace('.fits', '.noise.fits.gz'))

            # DC changing this for galaxy_plant()
            # segmap_file = "%s/%s" % (self.fake_image_path, fake_image_name.replace('.fits', '.check.fits'))
            segmap_file = "%s/%s" % (self.image_path, image_name.replace('.fits', '.check.fits'))

            #temp_field_match = re.findall(r"\d{2}_\d{3}[_]+\d{2}_\d{3}[_\d{1}\.]*", template_name)
            temp_field_match = template_name.split('.ut')[0]

            diff_dcmp_name = "{0}_{1}.diff.dcmp".format(fake_image_name[:-8],
                                                        template_name.replace(temp_field_match[:-1], "")[:-8])

            diff_dcmp_file = "%s/%s" % (self.diff_dir_path, diff_dcmp_name)

            filemissing = False
            for file in [image_file,
                         image_dcmp_file,
                         image_mask_file,
                         image_noise_file,
                         template_file,
                         template_dcmp_file]:

                if not os.path.exists(file):
                    filemissing = True
                    print('warning : file %s does not exist'%file)
            if filemissing: continue



            sexcat_good = None
            pixel_good = None
            if self.options.plant_in_galaxies:
                gal_demo_base_dir = "%s/%s" % (self.log_path, "galaxy_demographics")
                sexcat_files = glob.glob("%s/*_sexcat_good.txt" % gal_demo_base_dir)

                for s in sexcat_files:
                    sexcat_table = at.Table.read(s, format='ascii.ecsv')
                    model_props = sexcat_table.meta['comment']
                    f = model_props[0].split("=")[1].strip()


                    if f == image_file:
                        print("Found %s!" % f)
                        sexcat_good = s
                        break
            if sexcat_good is not None:
                pixel_good = sexcat_good.replace("sexcat", "pixel")


            file_associations[image_name] = FileAssociation(image_file,
                                                            image_dcmp_file,
                                                            image_mask_file,
                                                            image_noise_file,
                                                            fake_image_file,
                                                            fake_image_dcmp_file,
                                                            fake_image_mask_file,
                                                            fake_image_noise_file,
                                                            template_file,
                                                            template_dcmp_file,
                                                            diff_dcmp_file,
                                                            segmap_file,
                                                            sexcat_good,
                                                            pixel_good)

        return file_associations

    def append_fake_mag_file(self, line):
        with open(self.fake_mag_file, 'a') as fout:
            print(line, file=fout)

    def plant_fakes(self, fake_mag_range, clobber=False):

        # Get list of files to add fake stars to...
        psf_shape = 31


        # image_files = self.image_files
        # dcmp_files = self.get_matching_dcmp(image_files)
        # out_files = self.get_out_files(image_files)
        # mask_files = self.get_mask_files(image_files)

        # Plant fake stars
        for i, img in enumerate(self.image_names):

            print("Opening: %s" % img)
            try:
                file_association = self.file_associations[img]
            except:
                continue

            if not clobber and os.path.exists(file_association.fake_image_file):
                continue

            if self.options.plant_in_galaxies:
                print("running SExtractor")
                sextable = runsex(file_association.image_file, segmapname=file_association.segmap_file,
                                  zpt=fits.getval(file_association.image_file, 'ZPTMAG'))

                image_id = file_association.image_file.split('/')[-1].split('.')[0]
                glade_files = glob.glob('TileStats/*%s*txt'%image_id)
                for gf in glade_files:
                    #if len(glade_files):
                    #	glade_file = glade_files[0]
                    #else:
                    #	raise RuntimeError('temporarily raising error here! can\'t find GLADE file!')
                    good_ids, glade_ids, glade_bmags = np.array([]), np.array([]), np.array([])

                    try:
                        glade = at.Table.read(gf, format='ascii.ecsv')
                    except:
                        continue

                    for j in range(len(sextable.NUMBER)):
                        sep = astCoords.calcAngSepDeg(glade['Galaxy_RA'],
                                                      glade['Galaxy_Dec'],
                                                      sextable.X_WORLD[j],
                                                      sextable.Y_WORLD[j])*3600.
                        if not len(sep):
                            break
                        if min(sep) < 2:
                            good_ids = np.append(good_ids, sextable.NUMBER[j])
                            glade_ids = np.append(glade_ids, glade['Galaxy_ID'][sep == np.min(sep)][0])
                            glade_bmags = np.append(glade_bmags, glade['B'][sep == np.min(sep)][0])
                if not len(glade_ids):
                    continue

            image_hdu = fits.open(file_association.image_file)
            image_data = image_hdu[0].data.astype('float')

            # print("Opening: %s" % dcmp_files[i])
            dcmp_header = fits.getheader(file_association.image_dcmp_file)

            zpt = None
            try:
                zpt = dcmp_header['ZPTMAG']
            except KeyError:
                print("'ZPTMAG' keyword doesn't exist. Exiting...")
                continue
                #quit()

            psf_x, psf_xy, psf_y = dcmp_header['DPSIGX'],dcmp_header['DPSIGXY'], dcmp_header['DPSIGY']
            psf_model = self.generate_psf(psf_x, psf_xy, psf_y, psf_shape)

            fake_x = np.array([])
            fake_y = np.array([])
            fake_mags = np.random.uniform(fake_mag_range[0], fake_mag_range[1], fake_mag_range[2])
            gal_mags, glade_gal_mags, glade_gal_ids = [], [], []

            print("Opening: %s" % file_association.image_mask_file)
            mask_hdu = fits.open(file_association.image_mask_file)


            mask_data = mask_hdu[0].data.astype('float')
            if self.options.plant_in_galaxies:
                if not os.path.exists(file_association.segmap_file):
                    raise RuntimeError('segmap file %s was not created by SExtractor!' % file_association.segmap_file)

                segmap = fits.getdata(file_association.segmap_file)
                dcmp = txtobj(file_association.image_dcmp_file, cmpheader=True)

                # get rid of the stars
                self.options.glade = True

                if not self.options.glade:
                    for i in sextable.NUMBER[sextable.CLASS_STAR > 0.5]:
                        segmap[segmap == i] = 0
                else:
                    for i in sextable.NUMBER:
                        if i not in good_ids:
                            segmap[segmap == i] = 0

                for i in range(len(sextable.NUMBER)):

                    if not self.options.glade and sextable.CLASS_STAR[i] > 0.5:
                        continue
                    elif sextable.NUMBER[i] not in good_ids:
                        continue
                    if not self.options.glade:
                        sep = (dcmp.Xpos - sextable.X_IMAGE[i])**2. + (dcmp.Ypos - sextable.Y_IMAGE[i])**2.
                        sexmatch = np.where(sep == np.min(sep))[0]
                        if not len(sexmatch) or np.min(sep) > 5**2.:
                            segmap[segmap == sextable.NUMBER[i]] = 0
                            continue
                        sexmatch = sexmatch[0]

                        if dcmp.type[sexmatch] == '0x00000009' or dcmp.type[sexmatch] == '0x00000002' or dcmp.type[sexmatch] == '0x00000007':
                            segmap[segmap == sextable.NUMBER[i]] = 0

                good_indices = np.where((mask_data != 144.0) & (segmap != 0))
                #fits.writeto('segtmp.fits',segmap,overwrite=True)
                #import pdb; pdb.set_trace()
            else:
                good_indices = np.where(mask_data != 144.0)
            bad_indices = np.where(mask_data == 144.0)
            collisions = 0


            imshapex, imshapey = np.shape(image_data)
            #rand_ind_used = []
            #print('hi3')
            for j in range(len(fake_mags)):
                #print(j)
                candidate_found = False
                while not candidate_found:
                    if not len(good_indices[0]):
                        break

                    # Turns out that good_indices[1] corresponds to fits X-coord and good_indices[0] corresponds to fits Y-coord
                    rand_ind = random.choice(range(len(good_indices[1])))
                    #if rand_ind in rand_ind_used: continue
                    #else: rand_ind_used.append(rand_ind)
                    x = good_indices[1][rand_ind] + np.random.uniform(-0.5, 0.5)
                    y = good_indices[0][rand_ind] + np.random.uniform(-0.5, 0.5)

                    # Sanity -- Check against everything in the list for 10" separation
                    sep = (fake_x - x) ** 2. + (fake_y - y) ** 2.
                    #masksep = (bad_indices[1] - x)**2. + (bad_indices[0] - y)**2.

                    if x < 15 or y < 15 or x > imshapex-15 or y > imshapey - 15:
                        continue
                    #if min(masksep) < 30 ** 2: continue

                    if j == 0 or min(sep) > 5 ** 2:
                        fake_x = np.append(fake_x, x)
                        fake_y = np.append(fake_y, y)
                        if self.options.plant_in_galaxies:
                            try:
                                gal_mags.append(sextable.MAG_AUTO[sextable.NUMBER == segmap[int(y), int(x)]][0])
                            except:
                                gal_mags.append(-99)
                            if self.options.glade:
                                try:
                                    glade_gal_mags.append(glade_bmags[good_ids == segmap[int(np.round(y)),
                                                                                         int(np.round(x))]][0])
                                    glade_gal_ids.append(glade_ids[good_ids == segmap[int(np.round(y)),
                                                                                      int(np.round(x))]][0])
                                except: import pdb; pdb.set_trace()
                            segmap[int(y)-5:int(y)+6, int(x)-5:int(x)+6] = 0
                            good_indices = np.where((mask_data != 144.0) & (segmap != 0))
                        else:
                            gal_mags.append(-99)
                            glade_gal_mags.append(-99)
                            glade_gal_ids.append(-99)
                        candidate_found = True
                    else:
                        collisions += 1
                        # print("collision: (%s, %s); good fakes: %s" % (x,y, len(fake_x)))
            #print('hi4')
            print("# of collisions: %s" % collisions)
            # Build PSFs
            for m, g, gg, gi, x, y in zip(fake_mags, gal_mags, glade_gal_mags, glade_gal_ids, fake_x, fake_y):

                # Append log file with new fake point
                # gal_ra gal_dec psf_mag mag_gal mag_gal_GLADE GLADE_gal_id
                self.append_fake_mag_file('%s %s %s %s %s %s %s' % (file_association.image_dcmp_file, x, y, m, g, gg, gi))

                psf_mag = -2.5 * np.log10(np.sum(psf_model)) + zpt
                psf_flux = 10 ** (-0.4 * (m - psf_mag))

                max_size = np.shape(psf_model)[0]
                dx = dy = int((max_size - 1) / 2)

                image_data[int(y) - dy:int(y) + dy + 1, int(x) - dx:int(x) + dx + 1] \
                    += psf_model * psf_flux

            image_hdu[0].data[:] = image_data
            image_hdu.writeto(file_association.fake_image_file, clobber=True, output_verify='ignore')
            print('hi5')


    def galaxy_plant(self, gal_fake_mag_range, gal_fake_fwhm_factor, clobber=False):

        print("gal_fake_fwhm_factor: %s" % gal_fake_fwhm_factor)

        # Generate fakes for all files to be processed
        # fake_mags = np.random.uniform(gal_fake_mag_range[0], gal_fake_mag_range[1], gal_fake_mag_range[2])
        # fake_mags_iter = iter(fake_mags)

        # Loop over all files...
        for i, img in enumerate(self.image_names):
            print("Opening: %s" % img)
            try:
                file_association = self.file_associations[img]
            except:
                continue

            if not clobber and os.path.exists(file_association.fake_image_file):
                continue

            # Open the relevant data files, and extract vars
            image_hdu = fits.open(file_association.image_file)
            image_data = image_hdu[0].data.astype('float')
            segmap = fits.getdata(file_association.segmap_file)
            dcmp_header = fits.getheader(file_association.image_dcmp_file)
            zpt = None
            try:
                zpt = dcmp_header['ZPTMAG']
            except KeyError:
                print("'ZPTMAG' keyword doesn't exist. Exiting...")
                continue
            fwhm = dcmp_header['FWHM']
            pix_scale_arc_sec = np.abs(dcmp_header['CD2_2']) * 3600.0


            # Build the fake star model -- calculate PSF shape (e.g. the pixel-based radius for the fakes
            psf_x, psf_xy, psf_y = dcmp_header['DPSIGX'], dcmp_header['DPSIGXY'], dcmp_header['DPSIGY']

            psf_shape = int(np.ceil(gal_fake_fwhm_factor * fwhm))
            print("psf_shape: %s" % psf_shape)
            # must be odd...
            if psf_shape % 2 == 0:
                psf_shape += 1

            psf_model = self.generate_psf(psf_x, psf_xy, psf_y, psf_shape)
            psf_mag = -2.5 * np.log10(np.sum(psf_model)) + zpt
            dx = dy = int((psf_shape - 1) / 2)

            # We have already gone through the trouble of getting the "eligible" pixels. Rehydrate.
            pix_by_sexcat_id = {}
            galaxy_by_sexcat_id = {}

            # Keep track of the SExtractor/GLADE galaxy properties by sexcat_id
            sexcat_table = at.Table.read(file_association.sexcat_good, format='ascii.ecsv')
            sexcat_ids = list(sexcat_table["sexcat_id"])
            glade_ids = list(sexcat_table["glade_id"])
            glade_Bs = list(sexcat_table["glade_B"])
            sex_mags = list(sexcat_table["sex_mag"])
            for j, sid in enumerate(sexcat_ids):

                # Only process galaxies in our gal mag bin...
                if sex_mags[j] >= self.start_mag and sex_mags[j] <= self.end_mag:
                    galaxy_by_sexcat_id[sid] = (glade_ids[j], glade_Bs[j], sex_mags[j])

            # Keep track of the SExtractor galaxy pixels by sexcat_id
            pixel_table = at.Table.read(file_association.pixel_good, format='ascii.ecsv')
            sexcat_ids = list(pixel_table["sexcat_id"])
            good_pix_x = list(pixel_table["x"])
            good_pix_y = list(pixel_table["y"])
            for j, sid in enumerate(sexcat_ids):

                # Only process galaxies in our gal mag bin...
                if sid in galaxy_by_sexcat_id:
                    if sid not in pix_by_sexcat_id:
                        pix_by_sexcat_id[sid] = [[], []]

                    pix_by_sexcat_id[sid][0].append(good_pix_x[j])
                    pix_by_sexcat_id[sid][1].append(good_pix_y[j])


            # Actually inject the fakes
            injected_fakes = []
            for sid, good_pix in pix_by_sexcat_id.items():

                glade_id = galaxy_by_sexcat_id[sid][0]
                glade_B = galaxy_by_sexcat_id[sid][1]
                sex_mag = galaxy_by_sexcat_id[sid][2]

                all_x = good_pix[0]
                all_y = good_pix[1]

                x_min = np.min(all_x)
                x_max = np.max(all_x)
                y_min = np.min(all_y)
                y_max = np.max(all_y)

                psf_loc_x = np.arange(x_min, x_max, psf_shape)
                psf_loc_y = np.arange(y_min, y_max, psf_shape)

                for x in psf_loc_x:
                    for y in psf_loc_y:
                        if segmap[y, x] == sid:
                            # m = next(fake_mags_iter)
                            m = np.random.uniform(gal_fake_mag_range[0], gal_fake_mag_range[1])
                            injected_fakes.append((x, y, m, glade_id, glade_B, sex_mag))

            for x, y, m, gi, gg, g in injected_fakes:
                # import pdb; pdb.set_trace()

                psf_flux = 10 ** (-0.4 * (m - psf_mag))
                self.append_fake_mag_file('%s %s %s %s %s %s %s' % (file_association.image_dcmp_file, x, y, m, g, gg, gi))
                image_data[int(y) - dy:int(y) + dy + 1, int(x) - dx:int(x) + dx + 1] += psf_model * psf_flux

            image_hdu[0].data[:] = image_data
            image_hdu.writeto(file_association.fake_image_file, clobber=True, output_verify='ignore')

            # Write out Fakes for image file.
            fake_radius = (psf_shape / 2.0) * pix_scale_arc_sec
            fakes_region = "%s/%s" % (self.fake_log_path, img.replace("fits", "reg"))
            with open(fakes_region, 'w') as csvfile:

                csvfile.write("# Region file format: DS9 version 4.0 global\n")
                csvfile.write("global color=green\n")
                csvfile.write("image\n")

                for x, y, m, gi, gg, g in injected_fakes:
                    csvfile.write('circle(%s,%s,%s") # \n' % (x, y, fake_radius))
                    csvfile.write('point(%s,%s) # point=cross\n' % (x, y))

            print("Done w/ %s" % fakes_region)



    def generate_psf(self, psf_x, psf_xy, psf_y, psf_size):

        psf_model = np.zeros([psf_size, psf_size])
        for i in range(psf_size):
            for j in range(psf_size):
                x = i - psf_size/2
                y = j - psf_size/2
                zsq = 1/2.*(x**2./psf_x**2. + 2*psf_xy*x*y + y**2./psf_y**2.)
                psf_model[j, i] = (1 + zsq + 1/2.*zsq**2. + 1/6.*zsq**3.)**(-1)

        return psf_model

    def do_phot(self, iteration):

        # Copy over log files into the new fake log dir
        files = glob.glob('%s/*' % self.log_path)
        for f in files:
            if 'fake' not in f:
                os.system('cp %s %s' % (f, f.replace(self.image_dir, self.fake_image_dir)))

            # tap into this right after ABSPHOT, and before -diff MATCHTEMPL
            if 'ABSPHOT.outlist' in f:

                fin = open(f.replace(self.image_dir, self.fake_image_dir))
                fout = open(f.replace(self.image_dir, self.fake_image_dir).replace('.outlist', '.tmp.outlist'),'w')

                for line in fin:
                    line = line.replace('\n', '')

                    date_match = re.findall(r"ut1\d{5}", line)[0]

                    if self.image_dir == date_match:
                        lineparts = line.split()
                        lineparts[0] = lineparts[0].replace(date_match, '{0}_fake'.format(date_match))
                        lineparts[5] = lineparts[5].replace(self.image_dir, self.fake_image_dir)
                        lineparts[6] = lineparts[6].replace(self.image_dir, self.fake_image_dir)

                        print(" ".join(lineparts),file=fout)
                    else:
                        print(line.replace(date_match, '{0}_fake'.format(date_match)).replace(self.image_dir, self.fake_image_dir), file=fout)

                fout.close()
                fin.close()
                os.system('mv %s %s' % (f.replace(self.image_dir, self.fake_image_dir).replace('.outlist', '.tmp.outlist'),
                                        f.replace(self.image_dir, self.fake_image_dir)))

        for img in self.image_names:
            try: file_association = self.file_associations[img]
            except: continue
            os.system('rsync -avz %s %s' % (file_association.image_mask_file, file_association.fake_image_mask_file))
            os.system('rsync -avz %s %s' % (file_association.image_noise_file, file_association.fake_image_noise_file))
            os.system('rsync -avz %s %s' % (file_association.image_dcmp_file, file_association.fake_image_dcmp_file))

        os.system('pipeloop.pl -diff %s %s 1 -redo -stage MATCHTEMPL,DIFFIM,DIFFIMSTATS,DIFFDOPHOT,PIXCHK,DIFFCUT -k DC_MAX_NUMBER_OBJECTS 2000' % (self.fake_image_dir, self.template_dir))

    def get_phot(self, fake_mag_range, gal_mag_range):

        fake_mags = txtobj(self.fake_mag_file)
        for dcmp_path in np.unique(fake_mags.dcmpfile):
            # get key
            file_path = dcmp_path.replace('dcmp', 'fits')
            file_name = os.path.basename(file_path)
            try: file_association = self.file_associations[file_name]
            except: continue

            if not os.path.exists(file_association.diff_dcmp_file):
                print('Warning : diff dcmp file %s does not exist!' % file_association.diff_dcmp_file)
                badpart = file_association.diff_dcmp_file.split('.')[-3]
                dcmp_file_test = glob.glob(file_association.diff_dcmp_file.replace(badpart,'*'))
                if len(dcmp_file_test):
                    file_association.diff_dcmp_file = dcmp_file_test[0]
                else:
                    print('Warning : couldn\'t fix!  diff dcmp file %s does not exist!' % file_association.diff_dcmp_file)
                    continue

            dcmp = txtobj(file_association.diff_dcmp_file, cmpheader=True)
            zpt = fits.getval(file_association.diff_dcmp_file, 'IZPTMAG')
            for x, y, m, gm in zip(fake_mags.x[fake_mags.dcmpfile == file_association.image_dcmp_file],
                               fake_mags.y[fake_mags.dcmpfile == file_association.image_dcmp_file],
                               fake_mags.mag[fake_mags.dcmpfile == file_association.image_dcmp_file],
                               fake_mags.mag_gal[fake_mags.dcmpfile == file_association.image_dcmp_file]):

                sep = np.sqrt((x - dcmp.Xpos) ** 2. + (y - dcmp.Ypos) ** 2.)
                #if 's005aae22817' in file_association.diff_dcmp_file and abs(x-2816.96) < 0.02 and abs(y-824.21) < 0.01:
                #	import pdb; pdb.set_trace()
                    # and


                if len(np.where(sep < 2)[0]):
                    iMin = np.where(sep == np.min(sep))[0]
                    if len(iMin) > 1:
                        continue

                    # Sanity on the values...
                    if (dcmp.flux[iMin] < 0.0):
                        fout = open(self.out_match_file, 'a')
                        print('%s %.2f %.2f %.3f %.3f %.3f %.3f' % (file_association.diff_dcmp_file, x, y, m, gm, -99, -99), file=fout)
                        fout.close()
                    else:
                        fout = open(self.out_match_file, 'a')
                        print('%s %.2f %.2f %.3f %.3f %.3f %.3f' % (file_association.diff_dcmp_file, x, y, m, gm, -2.5 * np.log10(dcmp.flux[iMin]) + zpt, 1.086 * dcmp.dflux[iMin] / dcmp.flux[iMin]), file=fout)
                        fout.close()
                else:
                    fout = open(self.out_match_file, 'a')
                    print('%s %.2f %.2f %.3f %.3f %.3f %.3f' % (file_association.diff_dcmp_file, x, y, m, gm, -99, -99), file=fout)
                    fout.close()

            self.write_efficiency(fake_mag_range, file_name, file_association.diff_dcmp_file)
            if self.options.plant_in_galaxies:
                self.write_2d_efficiency(fake_mag_range, gal_mag_range, file_name, file_association.diff_dcmp_file)


        return self.out_match_file

    def get_phot_alliter(self, fake_mag_range, gal_mag_range, niter, out_match_files):

        fake_mags = txtobj(self.fake_mag_file)
        for dcmp_path in np.unique(fake_mags.dcmpfile):
            # get key
            file_path = dcmp_path.replace('dcmp', 'fits')
            file_name = os.path.basename(file_path)
            file_association = self.file_associations[file_name]

            if not os.path.exists(file_association.diff_dcmp_file):
                print('Warning : diff dcmp file %s does not exist!' % file_association.diff_dcmp_file)
                continue

            self.write_all_efficiency(fake_mag_range, file_name, file_association.diff_dcmp_file,out_match_files, niter)
            if self.options.plant_in_galaxies:
                self.write_all_2d_efficiency(fake_mag_range, gal_mag_range, file_name, file_association.diff_dcmp_file, out_match_files, niter)

    def write_efficiency(self, fake_mag_range, file_name, diff_dcmp_file):

        magbins=np.linspace(fake_mag_range[0], fake_mag_range[1], (fake_mag_range[1] - fake_mag_range[0])/fake_mag_range[3])

        # print(field_name)
        # print(diff_dcmp)

        om = txtobj(self.out_match_file)

        with open('%s/%s' % (self.efficiency_path, file_name.replace('.fits','.eff.txt')), 'w') as fout:
            print('# mag deteff N', file=fout)

            mag_sims = om.mag_sim[om.dcmpfile == diff_dcmp_file]
            detmags = om.detmag[om.dcmpfile == diff_dcmp_file]
            for m1, m2 in zip(magbins[:-1], magbins[1:]):

                iDet = np.where((mag_sims >= m1) & (mag_sims <= m2) & (detmags != -99))[0]
                iAll = np.where((mag_sims >= m1) & (mag_sims <= m2))[0]

                if len(iAll):
                    print('%.2f %.3f %s' % ((m1 + m2) / 2., len(iDet) / float(len(iAll)), len(iAll)), file=fout)
                else:
                    print('%.2f %.3f %s' % ((m1 + m2) / 2., np.nan, len(iAll)), file=fout)

    def write_2d_efficiency(self, fake_mag_range, gal_mag_range, file_name, diff_dcmp_file, doplot=True):

        galmagbins=np.linspace(gal_mag_range[0], gal_mag_range[1], (gal_mag_range[1] - gal_mag_range[0])/gal_mag_range[2])
        magbins=np.linspace(fake_mag_range[0], fake_mag_range[1], (fake_mag_range[1] - fake_mag_range[0])/fake_mag_range[3])

        # print(field_name)
        # print(diff_dcmp)

        om = txtobj(self.out_match_file)

        with open('%s/%s' % (self.efficiency_path, file_name.replace('.fits','.eff_gal.txt')), 'w') as fout:
            print('# mag gal_mag deteff N', file=fout)

            mag_sims = om.mag_sim[om.dcmpfile == diff_dcmp_file]
            mag_gal = om.mag_gal[om.dcmpfile == diff_dcmp_file]
            detmags = om.detmag[om.dcmpfile == diff_dcmp_file]
            for gm1, gm2 in zip(galmagbins[:-1], galmagbins[1:]):
                for m1, m2 in zip(magbins[:-1], magbins[1:]):

                    iDet = np.where((mag_sims >= m1) & (mag_sims <= m2) & (mag_gal >= gm1) & (mag_gal <= gm2) & (detmags != -99))[0]
                    iAll = np.where((mag_sims >= m1) & (mag_sims <= m2) & (mag_gal >= gm1) & (mag_gal <= gm2))[0]

                    if len(iAll):
                        print('%.2f %.2f %.3f %s' % ((m1 + m2) / 2., (gm1 + gm2) / 2., len(iDet) / float(len(iAll)), len(iAll)), file=fout)
                    else:
                        print('%.2f %.2f %.3f %s' % ((m1 + m2) / 2., (gm1 + gm2) / 2., np.nan, len(iAll)), file=fout)

        if doplot:
            eff = txtobj('%s/%s' % (self.efficiency_path, file_name.replace('.fits','.eff_gal.txt')))
            mag_grid = np.unique(eff.mag)
            gal_mag_grid = np.unique(eff.gal_mag)
            efficiency = np.zeros([len(mag_grid),len(gal_mag_grid)])
            for mg in range(len(mag_grid)):
                for gmg in range(len(gal_mag_grid)):
                    efficiency[mg,gmg] = eff.deteff[(eff.gal_mag == gal_mag_grid[gmg]) & (eff.mag == mag_grid[mg])]

            #efficiency = eff.deteff.reshape([len(mag_grid),len(gal_mag_grid)])
            import matplotlib.pyplot as plt
            plt.clf()
            ax = plt.axes()
            ax.imshow(efficiency,interpolation='nearest')
            ax.xaxis.set_ticks(range(len(gal_mag_grid)))
            ax.yaxis.set_ticks(range(len(mag_grid)))
            ax.xaxis.set_ticklabels(gal_mag_grid)
            ax.yaxis.set_ticklabels(mag_grid)
            ax.set_xlabel('Galaxy Mag')
            ax.set_ylabel('PSF Mag')

            plt.savefig('%s/%s' % (self.efficiency_path, file_name.replace('.fits','.eff.png')))
            plt.ion()
            plt.show()

    def write_all_efficiency(self, fake_mag_range, file_name, diff_dcmp_file, out_match_files, niter):
        from ast import literal_eval

        magbins=np.linspace(fake_mag_range[0], fake_mag_range[1], (fake_mag_range[1] - fake_mag_range[0])/fake_mag_range[3])

        om = txtobj(out_match_files[0])
        for omf in out_match_files[1:]:
            om_next = txtobj(omf)
            for k in om.__dict__.keys():
                om.__dict__[k] = np.append(om.__dict__[k],om_next.__dict__[k])

        with open('%s/%s' % (self.full_efficiency_path, file_name.replace('.fits','.eff.txt')), 'w') as fout:
            print('# mag deteff N', file=fout)

            galstr = 'gal_'

            iGood = "["
            for n in np.arange(niter)+1:
                diff_dcmp_iter = diff_dcmp_file.replace(self.fake_image_dir,"{0}_fake_{1}{2}".format(self.image_dir, galstr, n))
                iGood += "(om.dcmpfile == '%s') |"%diff_dcmp_iter
            iGood = eval(iGood[:-1]+"]")

            mag_sims = om.mag_sim[iGood] #[om.dcmpfile == diff_dcmp_file]
            detmags = om.detmag[iGood] #[om.dcmpfile == diff_dcmp_file]
            for m1, m2 in zip(magbins[:-1], magbins[1:]):

                iDet = np.where((mag_sims >= m1) & (mag_sims <= m2) & (detmags != -99))[0]
                iAll = np.where((mag_sims >= m1) & (mag_sims <= m2))[0]

                if len(iAll):
                    print('%.2f %.3f %s' % ((m1 + m2) / 2., len(iDet) / float(len(iAll)), len(iAll)), file=fout)
                else:
                    print('%.2f %.3f %s' % ((m1 + m2) / 2., np.nan, len(iAll)), file=fout)

    def write_all_2d_efficiency(self, fake_mag_range, gal_mag_range, file_name, diff_dcmp_file, out_match_files, niter, doplot=True):

        galmagbins=np.linspace(gal_mag_range[0], gal_mag_range[1], (gal_mag_range[1] - gal_mag_range[0])/gal_mag_range[2])
        magbins=np.linspace(fake_mag_range[0], fake_mag_range[1], (fake_mag_range[1] - fake_mag_range[0])/fake_mag_range[3])

        # print(field_name)
        # print(diff_dcmp)

        om = txtobj(out_match_files[0])

        for omf in out_match_files[1:]:
            om_next = txtobj(omf)
            for k in om.__dict__.keys():
                om.__dict__[k] = np.append(om.__dict__[k],om_next.__dict__[k])


        with open('%s/%s' % (self.full_efficiency_path, file_name.replace('.fits','.eff_gal.txt')), 'w') as fout:
            print('# mag gal_mag deteff N', file=fout)

            galstr = 'gal_'

            iGood = "["
            for n in np.arange(niter)+1:
                diff_dcmp_iter = diff_dcmp_file.replace(self.fake_image_dir,"{0}_fake_{1}{2}".format(self.image_dir, galstr, n))
                iGood += "(om.dcmpfile == '%s') |"%diff_dcmp_iter
            iGood = eval(iGood[:-1]+"]")

            mag_sims = om.mag_sim[iGood] #[om.dcmpfile == diff_dcmp_file]
            mag_gal = om.mag_gal[iGood] #[om.dcmpfile == diff_dcmp_file]
            detmags = om.detmag[iGood] #[om.dcmpfile == diff_dcmp_file]
            for gm1, gm2 in zip(galmagbins[:-1], galmagbins[1:]):
                for m1, m2 in zip(magbins[:-1], magbins[1:]):

                    iDet = np.where((mag_sims >= m1) & (mag_sims <= m2) & (mag_gal >= gm1) & (mag_gal <= gm2) & (detmags != -99))[0]
                    iAll = np.where((mag_sims >= m1) & (mag_sims <= m2) & (mag_gal >= gm1) & (mag_gal <= gm2))[0]

                    if len(iAll):
                        print('%.2f %.2f %.3f %s' % ((m1 + m2) / 2., (gm1 + gm2) / 2., len(iDet) / float(len(iAll)), len(iAll)), file=fout)
                    else:
                        print('%.2f %.2f %.3f %s' % ((m1 + m2) / 2., (gm1 + gm2) / 2., np.nan, len(iAll)), file=fout)

        if doplot:
            eff = txtobj('%s/%s' % (self.full_efficiency_path, file_name.replace('.fits','.eff_gal.txt')))
            mag_grid = np.unique(eff.mag)
            gal_mag_grid = np.unique(eff.gal_mag)
            efficiency = np.zeros([len(mag_grid),len(gal_mag_grid)])
            for mg in range(len(mag_grid)):
                for gmg in range(len(gal_mag_grid)):
                    efficiency[mg,gmg] = eff.deteff[(eff.gal_mag == gal_mag_grid[gmg]) & (eff.mag == mag_grid[mg])]

            #efficiency = eff.deteff.reshape([len(mag_grid),len(gal_mag_grid)])
            import matplotlib.pyplot as plt
            plt.clf()
            ax = plt.axes()
            ax.imshow(efficiency,interpolation='nearest',aspect=0.4)
            ax.xaxis.set_ticks(range(len(gal_mag_grid)))
            ax.yaxis.set_ticks(range(len(mag_grid)))
            ax.xaxis.set_ticklabels(gal_mag_grid)
            ax.yaxis.set_ticklabels(mag_grid)
            ax.set_xlabel('Galaxy Mag')
            ax.set_ylabel('PSF Mag')

            plt.savefig('%s/%s' % (self.efficiency_path, file_name.replace('.fits','.eff.png')))
            plt.ion()
            plt.show()
            print('%s/%s' % (self.efficiency_path, file_name.replace('.fits','.eff.png')))

    def hard_clean(self):
        self.soft_clean()

        lgdir = shutil.rmtree(self.fake_log_path.replace("/1", ""))
        try:
            print("'%s' Deleted..." % lgdir)
        except FileNotFoundError:
            print("'%s' not found. Proceeding..." % lgdir)

    def soft_clean(self):
        diff_log_dir = self.diff_dir_path.replace("workstch", "logstch")

        wkdir = self.diff_dir_path.replace("/1", "")
        try:
            shutil.rmtree(wkdir)
            print("'%s' Deleted..." % wkdir)
        except FileNotFoundError:
            print("'%s' not found. Proceeding..." % wkdir)

        stdir = self.fake_stitched_path.replace("/1", "")
        try:
            shutil.rmtree(stdir)
            print("'%s' Deleted..." % stdir)
        except FileNotFoundError:
            print("'%s' not found. Proceeding..." % stdir)

        diffwkdir = self.diff_dir_path.replace("/1", "")
        try:
            shutil.rmtree(diffwkdir)
            print("'%s' Deleted..." % diffwkdir)
        except FileNotFoundError:
            print("'%s' not found. Proceeding..." % diffwkdir)

        difflgdir = diff_log_dir.replace("/1", "")
        try:
            shutil.rmtree(difflgdir)
            print("%s Deleted..." % difflgdir)
        except FileNotFoundError:
            print("'%s' not found. Proceeding..." % difflgdir)

        try:
            print("'%s' Deleted..." % self.efficiency_path)
        except FileNotFoundError:
            print("'%s' not found. Proceeding..." % self.efficiency_path)

        try:
            os.remove(self.fake_mag_file)
            print("'%s' Deleted..." % self.fake_mag_file)
        except FileNotFoundError:
            print("'%s' not found. Proceeding..." % self.fake_mag_file)


class AllStages():
    def __init__(self):
        pass

    def add_options(self, parser=None, usage=None):

        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('--stage', default='all', type='string', help='Comma-separated stages. Options are plant, photpipe, getphot')
        parser.add_option('--clobber', default=False, action="store_true", help='overwrite previous results if set')

        parser.add_option('--root_path', default='$PIPE_DATA/workstch', type='string', help='Root directory')
        parser.add_option('--image_dir', default='gw190425', type='string', help='Image directory')
        parser.add_option('--field_name_start', default='s005', type='string', help='Image directory')
        parser.add_option('--template_dir', default='gw190425tmpl', type='string', help='Template directory')
        parser.add_option('--image_list', default='gw190425/<change this to the one you want>.txt', type='string', help='File with all observations')
        parser.add_option('--template_list', default='gw190425/<change this to the one you want>.txt', type='string', help='File with all templates')

        parser.add_option('--iterations', default=1, type='int', help='Number of times to run the job')
        parser.add_option('--iteration_start', default=1, type='int', help='Integer used in directory name to start the iterations')

        parser.add_option('--gal_bin_to_process', default='', type='string', help='Gal mag range for galaxy injections')
        parser.add_option('--gal_mag_range', default=(13, 22, 0.5), nargs=2, type='float', help='gal mag tuple: (min, max, bin size)')
        parser.add_option('--fake_mag_range', default=(18, 25, 1500, 0.2), nargs=2, type='float', help='Fake mag tuple: (min, max, # of stars, bin size)')
        parser.add_option('--gal_fake_mag_range', default=(18, 25, 5000, 0.2), nargs=2, type='float',help='Fake mag tuple: (min, max, # of stars, bin size)')
        parser.add_option('--gal_fake_fwhm_factor', default=2.0, type='float', help='Factor times image FWHM to set fake star psf size in pixels')

        parser.add_option('--plant_in_galaxies', default=False, action="store_true",
                          help='if set, generate positions from the SExtractor galaxy mask')
        parser.add_option('--clean', default=False, action="store_true",
                          help='if set, run hard clean')

        return(parser)


if __name__ == "__main__":

    # python DetEff_StaticFile.py --iterations 2 --image_dir gw190814 --template_dir gw190814tmpl --field_name_start s --image_list gw190814/all_tiles_ascii.txt --template_list gw190814/all_temps.txt

    # python DetEff_StaticFile.py --stage photpipe
    # python DetEff_StaticFile.py --stage plant --plant_in_galaxies --iterations 6
    # python DetEff_StaticFile.py --stage photpipe --image_dir gw170814_2 --template_dir gw170814tmpl --image_list all_tiles_ascii.txt --template_list all_temps.txt

    import optparse

    usagestring='USAGE: DetEff.py'

    allstages = AllStages()
    parser = allstages.add_options(usage=usagestring)
    options,  args = parser.parse_args()

    detEff = DetermineEfficiencies(root_path=os.path.expandvars(options.root_path),
                                    image_dir=options.image_dir,
                                    template_dir=options.template_dir,
                                    image_list=options.image_list,
                                    template_list=options.template_list)
    detEff.options = options

    iterations = options.iterations

    t0 = time.time()

    i = options.iteration_start
    out_match_files = []
    while i <= iterations:

        detEff.initialize(i, options.gal_bin_to_process)
        if options.clean:
            detEff.hard_clean()
            continue

        if 'plant' in options.stage or 'all' in options.stage:

            if options.plant_in_galaxies:
                if options.gal_bin_to_process == "":
                    raise Exception("Which gal mag bin to process?")
                detEff.galaxy_plant(gal_fake_mag_range=options.gal_fake_mag_range, gal_fake_fwhm_factor=options.gal_fake_fwhm_factor)
            else:
                detEff.plant_fakes(fake_mag_range=options.fake_mag_range, clobber=options.clobber)

        if 'photpipe' in options.stage or 'all' in options.stage:
            detEff.do_phot(iteration=i)

        if 'getphot' in options.stage or 'all' in options.stage:
            out_match_files += [detEff.get_phot(fake_mag_range=options.fake_mag_range, gal_mag_range=options.gal_mag_range)]

        i += 1

    if 'getphot' in options.stage or 'all' in options.stage:
        detEff.get_phot_alliter(fake_mag_range=options.fake_mag_range,gal_mag_range=options.gal_mag_range,niter=iterations,out_match_files=out_match_files)


    t1 = time.time()
    total = t1 - t0

    print("\n*****************\nProcess time: %s\n*****************\n" % total)
