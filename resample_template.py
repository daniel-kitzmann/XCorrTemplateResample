
import numpy as np
import matplotlib.pyplot as plt
import argparse
from astropy.io import fits
from os.path import exists
import sys



def parseCommandLineArgs():

  args_ok = True

  parser = argparse.ArgumentParser(
    description= "This script supports the following command-line parameters")

  parser.add_argument(
    'filename', 
    metavar='template file path',
    type=str, 
    nargs=1,
    help='file path for the template fits file')
  parser.add_argument(
    '-r', 
    dest='resolution', 
    help='spectral resolution', 
    required=True)
  parser.add_argument(
    '-o', 
    help='output file name',
    required=False)
  parser.add_argument(
    '-b', 
    help='wavelength bin in nm', 
    required=False)

  args = parser.parse_args()


  template_file_path = args.filename[0]

  resolution = float(args.resolution)


  min_wavelength = -1.0
  max_wavelength = -1.0

  if args.b:

    bin_boundaries = args.b.split(',')
  
    if len(bin_boundaries) != 2:
      print('The script requires two wavelength boundary points, separated by a comma!')
      args_ok = False

    bound_a = float(bin_boundaries[0])
    bound_b = float(bin_boundaries[1])
  
    if bound_a == bound_b:
      print('The wavelength boundary points can\'t be identical!')
      args_ok = False

    if bound_a < bound_b:
      min_wavelength = bound_a
      max_wavelength = bound_b
    else:
      min_wavelength = bound_b
      max_wavelength = bound_a

  if args.o:
    output_file_path = args.o
  else:
    file_path_comp = template_file_path.split('.fits')

    output_file_path = file_path_comp[0] + '_r_' + str(np.round(resolution)).split('.')[0]

    if len(file_path_comp) > 1:
      output_file_path = output_file_path + '.fits'

  return args_ok, template_file_path, output_file_path, min_wavelength, max_wavelength, resolution


#Reads the fits file and extracts the data
def readTemplateFile(file_path):
  
  template_wavelengths = np.array([])
  template = np.array([])
  template_max_wavelength = 0
  template_min_wavelength = 0

  if exists(file_path):
    print('Reading template file', file_path)

    fits_data = fits.open(file_path)
    template_wavelengths = np.array(fits_data[0].data[0,:])
    template =  np.array(fits_data[0].data[1,:])

    template_max_wavelength = template_wavelengths[0]
    template_min_wavelength = template_wavelengths[-1]
  else:
    print('Template file', file_path, 'does not exist!')

  return template_wavelengths, template, template_max_wavelength, template_min_wavelength



def saveFITSTemplate(file_path, template_wavelengths, template):
  
  print('Saving new template to:', file_path)

  data = np.stack((template_wavelengths, template))

  fits.writeto(file_path, data, overwrite=True)

  if exists(file_path) == False:
    print('Couldn\'t save output file to', file_path)

  return




def setBinBoundaries(wavelengths, template_wavelengths):
  
  bin_boundaries = np.zeros(wavelengths.size+1)

  bin_boundaries[0] = wavelengths[0] - (wavelengths[1] - wavelengths[0])*0.5
  bin_boundaries[-1] = wavelengths[-1] + (wavelengths[-1] - wavelengths[-2])*0.5

  for i in range (1,wavelengths.size):
    bin_boundaries[i] = wavelengths[i-1] + (wavelengths[i] - wavelengths[i-1])*0.5

  if bin_boundaries[0] > template_wavelengths[0]:
    bin_boundaries[0] = template_wavelengths[0]

  if bin_boundaries[-1] < template_wavelengths[-1]:
    bin_boundaries[-1] = template_wavelengths[-1]

  return bin_boundaries



#find the closest index in the tabulated template wavelength array
#this function assumes a constant wavenumber step
def findIndex(wavelength, template_wavelengths, side='left', wavenumber_step = 0.01):

  #starting wavenumber in cm-1 assuming that the wavelength is in nm
  nu_start = np.round(1.0/template_wavelengths[0] * 1e7, 6)
  
  #wavenumber of the wavelength we're seeking the index for
  nu = np.round(1.0/wavelength * 1e7, 6)

  #calculate the array index, assuming a constant wavenumber step
  index = int((nu - nu_start)/wavenumber_step)

  index_ret = index

  if side=='right':
    if template_wavelengths[index] < wavelength:
      index_ret = index
    else:
      index_ret = index+1
  else:
    if template_wavelengths[index] > wavelength:
      index_ret = index
    else:
      index_ret = index-1

  return index_ret



#create the new wavelength grid based on the desired resolution
def createWavelengthGrid(resolution, min_wavelength, max_wavelength, template_wavelengths):

  #first, we count the number of points we need for the pre-allocation
  wavelength = max_wavelength
  nb_points = 1

  while wavelength > min_wavelength:
    wavelength -= wavelength/resolution
    
    if wavelength < min_wavelength:
      break

    nb_points += 1


  #now we create the new wavelength grid with a constant resolution
  wavelengths = np.zeros(nb_points)
  wavelengths[0] = max_wavelength

  for i in range(1,nb_points):
    wavelengths[i] = wavelengths[i-1] - wavelengths[i-1]/resolution

  bin_boundaries = setBinBoundaries(wavelengths, template_wavelengths)


  #check the minumum resolution
  index = findIndex(max_wavelength, template_wavelengths, 'right')
  min_resolution = template_wavelengths[index] / (template_wavelengths[index] 
    - template_wavelengths[index+1])

  if (min_resolution < resolution):
    print('Warning: Desired resolution of', resolution, 
          'is higher than available in the template:',min_resolution, 
          'at the largest wavelength', max_wavelength, 'nm', 
          'This will lead to an oversampling of the original template data.')
  

  return wavelengths, bin_boundaries




def integrateBin(bin, center_wavelength, template_wavelengths, template):

  #look up the wavelength indices of the bin boundaries
  indices = np.array([findIndex(bin[0], template_wavelengths, 'right'), 
                      findIndex(bin[1], template_wavelengths, 'left')])

  #in case of oversampling, we go back to linear interpolation between existing data points
  if indices[0] == indices[1]:
     #find the index of closest data point to the bin center in the template
     index = findIndex(center_wavelength, template_wavelengths)

     #linear interpolation between two adjacent template data points
     bin_value = template[index] + (template[index+1] - template[index])/(
       (template_wavelengths[index+1] - template_wavelengths[index]))*(
       (center_wavelength - template_wavelengths[index]))
  else:
    #otherwise, we do a trapezoidal quadrature over the spectral bin to calculate the mean template value

    #linear interpolation of the bin boundaries
    left_boundary = template[indices[0]-1] + (template[indices[0]] - template[indices[0]-1]) \
      /(template_wavelengths[indices[0]] - template_wavelengths[indices[0]-1]) \
      * (bin[0] - template_wavelengths[indices[0]-1])
    right_boundary = template[indices[-1]] + (template[indices[-1]+1] - template[indices[-1]]) \
      /(template_wavelengths[indices[-1]+1] - template_wavelengths[indices[-1]]) \
      * (bin[1] - template_wavelengths[indices[-1]])

    bin_wavelengths = np.concatenate( (
      np.array([bin[0]]), 
      template_wavelengths[indices[0]:indices[1]], 
      np.array([bin[1]])))
    bin_template = np.concatenate(( 
      np.array([left_boundary]), 
      template[indices[0]:indices[1]], 
      np.array([right_boundary])))

    bin_value = np.trapz(bin_template, bin_wavelengths)/(bin[1] - bin[0])

  return bin_value




def main():

  args_ok, template_file_path, output_file_path, min_wavelength, max_wavelength, resolution = parseCommandLineArgs()

  if args_ok == False:
    exit()


  template_wavelengths, template, template_max_wavelength, template_min_wavelength = readTemplateFile(
    template_file_path)

  if template_wavelengths.size == 0:
    exit()


  #if no command line argument for re-binning is given, use the entire wavelength range
  if max_wavelength == -1 and min_wavelength == -1:
    max_wavelength = template_max_wavelength
    min_wavelength = template_min_wavelength

  #check the wavelength ranges
  if max_wavelength > template_max_wavelength or min_wavelength < template_min_wavelength:
    print('The supplied wavelength boundaries', max_wavelength, 
          'and', min_wavelength, 
          'are outside the template wavelength interval:', 
          template_max_wavelength, template_min_wavelength)
    exit()


  wavelengths, wavelength_bins = createWavelengthGrid(
    resolution, 
    min_wavelength, 
    max_wavelength, 
    template_wavelengths)

  #create the new template
  new_template = np.zeros(wavelengths.size)

  for i in range(wavelengths.size):
    percentage = float(i/wavelengths.size) * 100
    sys.stdout.write("Processing template: {:.2f}%\r".format(percentage))
    sys.stdout.flush()

    new_template[i] = integrateBin(
      np.array((
        wavelength_bins[i], 
        wavelength_bins[i+1])), 
        wavelengths[i], 
        template_wavelengths, 
        template)


  saveFITSTemplate(output_file_path, wavelengths, new_template)


  plt.plot(template_wavelengths, template)
  plt.plot(wavelengths, new_template)
  plt.show()



if __name__ == "__main__":
    main()




