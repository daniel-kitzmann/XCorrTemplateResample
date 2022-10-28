# Cross-Correlation Template Resampling #
#### Author: Daniel Kitzmann ####

The python script *resample_template.py* allows to resample the single-species
templates in their FITS format to a constant resolution and optionally also to 
a specific wavelength range. 

The script is called by:

python resample_template.py template_file.fits -r xxxxx

where *template_file.fits* is the path to the original template and *xxxxx* is 
the desired resolution. By default, the resampled template will be saved in the
file *template_file_r_xxxxx.fits* within the same file path as the original template.

The resampling is done by integrating the original high-resolution to a lower 
resolution. It is also possible to resample a template to a higher resolution 
than that provided by the original template. In this case, the script will 
resample using linear interpolation.

The script has also two additional, optional parameters:

-o filename

specifies a path where the resampled template should be saved

-b x y

specifies a limited wavelength range between *x* and *y* the templated should be 
resampled in. The parameters *x* and *y* are the wavelengths in units of nm.

