# USAGE: vmd -dispdev text -e generate_psf-pdb.tcl

# Load System PSF and PDB
mol new raw/step3_input.psf
mol addfile raw/step3_input.pdb

# Select Residues
set sel [atomselect top "resname ACN or name FE3P"]

# Write PSF and PDB
$sel writepsf aerobactin_fe.psf
$sel writepdb aerobactin_fe.pdb

exit
