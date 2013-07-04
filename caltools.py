# caltools.py

# A module is a Python object with 
# arbitrarily named attributes that you can bind and reference.

# to use: 
# import caltools

import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from scipy.constants import pi as pi
import math

'''
ROLO bands with "bad seeing" aka bad data... 
this is true for the images, but are the chips also flawed?
4,6,9,13,16,19
19,415 # bad seeing
4,442 # bad seeing
6,488 # bad seeing
16,347 # bad seeing
'''

def read_rolo_chips(root):
    '''
    Input: a path containing the file.
    Returns: a pd.DataFrame
    JULIAN_DATE, FILTER_ID, PHASE, SUB_EARTH_LAT, SUB_EARTH_LON, SUB_SOLAR_LAT, 
    SUB_SOLAR_LON, CALIBRATION, EXTINCTION, CHIPS
    where CHIPS are 11 entries that are pixel area averages in instrument units, 
    DN/sec. Names and orders of chips: MS2, E, N, S, W, Aris3, Aris7, Coper, 
    Tycho, Highl, TyRay

    Reference: Kieffer and Stone (2005) The Spectral Irradiance of the Moon,
    The Astronomical Journal, 129:2887-2901.

    The values for PHASE are geocentric lunar phase.

    Some of the chip values are negative, and some are valued -7.7700000E+02 
    -- this is an indicator that the chip area is either not illuminated or not
    visible. Also, there are some anomalous EXTINCTION values to look out for.
    '''

    columns = ['julian_date', 
                'filter_id', 
                'phase', 
                'subearth_lat', 
                'subearth_lon', 
                'subsolar_lat',
                'sub_solar_lon',
                'calibration',
                'extinction',
                'serenitatis',
                'east',
                'north',
                'south',
                'west',
                'aristarchus3',
                'aristarchus7',
                'copernicus',
                'tycho',
                'highlands',
                'tychoray'
    ]

    chips = pd.read_fwf(
        root + 'alex_311.chips32.asc', 
        widths   = (16, 5) + 5*(11,) + 13*(17,),
        skiprows = 0, 
        header   = None,
        names    = columns
    )

    return chips


def get_filter_lookup(path): # should go in toolkit
    '''
    Read in the FILTER_ID lookup table
    Return pd.DataFrame indexed to filter_id
    '''
    lookup = pd.read_csv(
        path + 'lookup.csv', 
        header = 0,
        index_col = 'filter_id' 
        )

    return lookup


def read_phase(root):
    '''
    Topocentric phase angles for the ROLO chip locations are in 
    'alex_311_phased32.asc', indexed to the chips data file 
    'alex_311.chips32.asc'. Use the topocentric phase angles
    Negative phase angle values correspond to waxing lunar phases.
    Is PHASE in degrees? SOLAR_DISTANCE is in AU
    Anomalous values correspond to times where the area was not illuminated
    '''

    phase = pd.read_fwf(
        root + 'alex_311_phased32.asc', 
        widths   = (16, 5, 12, 12), 
        header   = None, 
        names    = ('JULIAN_DATE', 'FILTER_ID', 'PHASE', 'SOLAR_DISTANCE'), 
        skiprows = 0
        )

    return phase


def read_angles(root):
    '''
    'angles_new_3.dat' contains incidnece, emission, and phase 
    in units of radians for the 11 ROLO chips

    Some of the angles are greater than 360 degrees, but modulo 360.0 gives 
    reasonable values.
    '''

    emission_col = ['serenitatis_emission', 'east_emission', 'north_emission',
                    'south_emission', 'west_emission', 'aristarchus3_emission', 
                    'aristarchus7_emission', 'copernicus_emission', 
                    'tycho_emission', 'highlands_emission', 'tychoray_emission'
                    ]
    phase_col = ['serenitatis_phase', 'east_phase', 'north_phase', 'south_phase',
                'west_phase', 'aristarchus3_phase', 'aristarchus7_phase',
                'copernicus_phase', 'tycho_phase', 'highlands_phase', 
                'tychoray_phase'
                ]
    incidence_col = ['serenitatis_incidence', 'east_incidence', 'north_incidence',
                    'south_incidence', 'west_incidence', 'aristarchus3_incidence',
                    'aristarchus7_incidence', 'copernicus_incidence',
                    'tycho_incidence', 'highlands_incidence', 'tychoray_incidence'
                    ]
    columns = incidence_col + emission_col + phase_col

    angles=pd.read_fwf(
        root + 'angles_new_3.dat',
        widths   = (12,) + 32*(13,),
        header   = None,
        skiprows = 0,
        names    = columns
        )
    angles = angles*(180/pi) # change radians to degrees
    angles = np.fmod(angles, 360) # TODO: I don't think this is working
    return angles


def interpolate_solar_spectrum(path, lookup):
    #TODO needs error catching for out-of-bounds input 
    '''
    Input: path and a pd.df with filter_id and band info

    Returns: Returns interpolated solar spectrum an np array
    Interpolates solar spectrum to a set of wavelengths bands [nm]
    
    Solar spectrum is in units of irradiance: 'W m^-2 nm^-1'
    '''
    # read in the SOLAR IRRADIANCE data file: 350 to 1100 nm
    solar = pd.read_csv(path+'solarspectrum.csv', header=1)
    # Interplolate the SOLAR IRRADIANCE to fit bands
    x,y = solar.wavelength, solar.solar_irradiance
    interpolation = interp1d(x, y, kind='cubic')
    # filter_id, band, solar_irradiance
    bands = lookup.sort(['band_nm'])['band_nm'].values # array of bands
    interp_solar = interpolation(bands[1:24])
    interp_solar_df = pd.DataFrame(interp_solar, 
        index=bands[1:24], columns=['solar_irradiance'])

    return interp_solar_df


def get_chip_reflectance(path, chips, lookup, feature): 
    '''
    Input: 
    Takes ROLO chip data (1 chip/roi) and returns I/F (reflectance)
    IOF = (pi * RADIANCE) / INTERPOLATED_SOLAR_IRRADIANCE
    Looks at the filter_id and then multiplies by the appropriate
    solar irradiance.
    '''

    # interpolate solar spectrum to rolo bands (returns pd.df)
    interp_solar_df = interpolate_solar_spectrum(path, lookup)

    bands = lookup.sort(['band_nm'])['band_nm'].iloc[1:24] #subset of lookup

    # You can do math on series... but if a selection results in a new pd.df, 
    # doing math on that won't change the values in the original df
    # you also can't assign a varaible to the series
    for index in bands.index:
        # chips[feature][chips['filter_id'] == index] # a series
        # interp_solar is indexed by band
        band = int(lookup.loc[index].values)
        # What happens when 'band' is not in interp_solar_df? KeyError
        chips[feature][chips['filter_id'] == index] = (chips[feature][chips['filter_id'] == index]*pi)/interp_solar_df.loc[band].values
    
    return chips


def interp2wacbands(bands, data):
    '''
    Input: bands and data for what you want to interpolate
    Returns: data interpolated to LROC WAC bands
    '''
    wac_bands = [360, 415, 566, 604, 643, 689] # removed 320 nm band
    x,y = bands, data
    interpolation = interp1d(x, y, kind='cubic')
    return interpolation(wac_bands)


def lommel_seeliger(incidence, emission):
    '''
    Calculate Lommel-Seeliger function
    Input: emission (e) and incidence (i) angles
    Lommel-Seeliger = cos(i) / (cos(e) + cos(i)) 

    From Hapke et al., 2012:
    "Because theoretical photometric functions of particulate media based on 
    the equation of radiative transfer contain the Lommel-Seeliger function 
    as a common factor, all values of I/F were divided by LS."

    Citation: Fairbairn, M. B. (2005) "Planetary Photometry: The Lommel-Seeliger 
    Law" Journal of the Royal Astronomical Society of Canada, p. 92-93.
    '''

    cos_incidence = np.cos(np.radians(incidence))
    cos_emission = np.cos(np.radians(emission)) 
    LS = cos_incidence / (cos_emission + cos_incidence)
    return LS


'''
Cool trick:

b=( a*(b^2)) / (abs(ab))

if a is negative and b is positive, b becomes negative
if a is negative and b is negative, then b stays negative
if a is positive and b is positive then nothing happens
'''