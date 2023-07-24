
from astropy.io import fits
from astropy.table import Table

on_att_filename = 'analysis/attitude.fits'
hdu_att_list = fits.open(on_att_filename)
att_entries = Table.read(on_att_filename, hdu=1)
att_ra  = att_entries[int(len(att_entries)/2)]['AHFRA']
att_dec = att_entries[int(len(att_entries)/2)]['AHFDEC']

outfile_path = 'analysis/pointing_coord.txt'
with open(outfile_path, 'w') as file:
    file.write('circle %s %s 1.0\n'%(att_ra,att_dec))
    file.write('circle %s %s 1.0\n'%(att_ra,att_dec+0.2))
    file.write('circle %s %s 1.0\n'%(att_ra+0.2,att_dec))
    file.write('circle %s %s 1.0\n'%(att_ra-0.2,att_dec))
    file.write('circle %s %s 1.0\n'%(att_ra,att_dec-0.2))

