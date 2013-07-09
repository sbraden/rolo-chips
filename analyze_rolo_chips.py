# python
# Sarah Braden
# 7 May 2013

import pandas as pd
import numpy as np
import math
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from scipy.constants import pi as pi
import itertools
import caltools # my own module


colors = [
  'red',
  'blue',
  'green',
  'burlywood',
  'cyan',
  'magenta',
  'yellow',
  'black',
  'darkorange',
  'purple',
  'lightgreen'
]

colorloop=intertools.cycle(colors) # using intertools!


def read_WAC(root, feature): # wac_bands used elsewhere too
    '''
    Read in WAC data with 2 files: image_data and spectra_data
    'spectra_data' is in units of calibrated DN.
    '''
    spectra = pd.read_csv(
        root + feature +'_WAC_spectra.csv', 
        header   = 0,
        index_col = 0
    )
    image_data = pd.read_csv(
        root + feature +'_WAC_image_data.csv', 
        header   = 0,
        index_col = 0
    )
    # add column from one df to another (with the same index)
    spectra['exp_time'] = pd.Series(image_data['exp_time'].values, 
                                    index=spectra.index
    )
    # add column from one df to another (with the same index)
    spectra['subsolar_ground_azimuth'] = pd.Series(
                            image_data['subsolar_ground_azimuth'].values, 
                            index=spectra.index
    )

    wac_bands = [360, 415, 566, 604, 643, 689]
    columns = ['band','phase','avg','stdev','num'] # this is the same as before
    phase_groups = [15, 20, 25, 30, 35, 40] # same as the other phase_groups array

    subspectra1 = spectra[spectra.subsolar_ground_azimuth < 180]
    subspectra2 = spectra[spectra.subsolar_ground_azimuth > 180]
    x=subspectra1['avg_pha'].values
    y=subspectra1['689'].values/subspectra1['exp_time'].values
    plt.plot(x, y, 'or')
    x=subspectra2['avg_pha'].values
    y=subspectra2['689'].values/subspectra2['exp_time'].values
    plt.plot(x, y, 'ob')
    plt.xlabel('phase [degrees]', fontsize=14)
    plt.ylabel('DN/ms', fontsize=14)
    plt.title('WAC_serenitatis blue: >180, red: <180 azimuth')
    plt.show()
    #y=spectra['avg_pha'].values


    data = []
    for phase in phase_groups:
        group = spectra[(spectra['avg_pha'] > phase-2) & (spectra['avg_pha'] < phase+2)]
        for band in wac_bands:
            cal_dn = group[str(band)].values/group['exp_time'].values
            output = band, phase, cal_dn.mean(), cal_dn.std(), len(cal_dn)
            data.append(output)
    cal_dn_df = pd.DataFrame(data, columns=columns)
    grouped = cal_dn_df.groupby(['phase', 'band'])

    return grouped.aggregate(np.mean) # Giving correct values?


def main():

    WAC_root = '/home/sbraden/data/rolo_chip_cal/'

    ROLO_root = '/home/sbraden/Dropbox/calibration/ROLO/datasetChips/'
    feature = 'serenitatis' 

    # Read in ROLO data and process into I/F reflectance
    lookup = caltools.get_filter_lookup(roloroot) # filter_id => band lookup table
    chips = caltools.process_rolo_data(roloroot, feature, lookup)

    # STEP 3: (optional) appy the Lommel-Seeliger correction
    incidence = chips[feature+'_incidence']
    emission = chips[feature+'_emission']
    LS = caltools.lommel_seeliger(incidence, emission)
    chips['LS'] = pd.Series(LS.values)
    #chips[feature] = chips[feature]*LS # apply LS

    # STEP 4: set data to only the appropriate phases by setting others to nan
    # Negative phase angle values correspond to waxing lunar phases
    chips.topo_phase[chips.topo_phase > 0] = np.nan
    chips.topo_phase = abs(chips.topo_phase.values)

    ###########################################

    # groups for phase angles
    phase_groups = [15, 20, 25, 30, 35, 40]
    # STEP 3: all values for the chip of interest are in I/F, 
    # index chips by filter_id
    # df=pd.DataFrame(chips, index=chips['filter_id'])
    # sort by phase, bin by phase?, then group/aggregate, interp to WAC bands

    # see if I can replace this loop by some other kind of selection
    data = []
    columns = ['band','phase','avg','stdev','num']

    for phase in phase_groups:
        group = subchip[(subchip['real_phase'] > phase-2) & (subchip['real_phase'] < phase+2)]
        output = band, phase, chip_rad.mean(), chip_rad.std(), len(chip_rad.index)
        data.append(output)

    all_data = pd.DataFrame(data, columns=columns)
    grouped = all_data.groupby(['phase', 'band']) # indexed by phase group then band
    # print all_data.groupby(['band', 'phase']).mean().avg
    # print all_data.groupby(['band', 'phase']).mean().std
 
    # create a MultiIndex from the groups through Aggregation
    # This averages multiple observations in the same filter together
    avg_multi_index = grouped.aggregate(np.mean)
    averages = avg_multi_index['avg']
    stdevs = avg_multi_index['stdev']

    # READ IN AND PROCESS WAC Data
    wac_dataframe = read_WAC(WAC_root, feature)
    # print wac_dataframe

    # not the best way to define this... ok for now
    #customized
    rolo_bands = [353, 405, 413, 467, 476, 550, 545, 555, 667, 695, 706]

    for phase in phase_groups:
        # interpolate data and stdev to WAC bands
        rolo_iof_interp = caltools.interp2wacbands(rolo_bands, iof_df.values)
        rolo_iof_stdev_interp = caltools.interp2wacbands(rolo_bands,
            iof_df_stdev.values)
        cal_factor = wac_dataframe.loc[phase].avg/rolo_iof_interp

        x = cal_factor.index
        y = cal_factor.values
        p = rolo_iof_interp
        p_stdev = rolo_iof_stdev_interp
        r = wac_dataframe.loc[phase].avg
        r_stdev = wac_dataframe.loc[phase].stdev
        u = np.power(p_stdev/p, 2)
        w = np.power(r_stdev/r, 2)
        yerr = y * np.power(u+w, 0.5)
        plt.errorbar(x, y, yerr=yerr, marker='o', label=str(phase), 
            c=colorloop.next())
        # print phase, cal_factor

    plt.legend(loc='best')
    plt.xlabel('WAC bands [nm]', fontsize=14)
    plt.ylabel('ROLO-derived calibration factor', fontsize=14)
    filename = feature+ '_cal_factor_abs'
    plt.title(filename)
    plt.savefig(filename+'.png', dpi=200)
    plt.close

    # test_multi_index.index.names # gives names of indices

    #plot_and_save(
    #    band_chip[2], 
    #    chip_radiance, 
    #    'geocentric phase', 
    #    'ROLO 550 nm radiance', 
    #    '550nm_geocentricphasetest')

if __name__ == '__main__':
  main()

'''
JUNK

    # match input to solar irrandiance data so interp won't throw error
    solarmax = solar['wavelength'].values.max()
    solarmin = solar['wavelength'].values.min()
    x = np.where((bands>solarmin) & (bands<solarmax)) # returns indices
    bands4interp = bands[x[0]] # return values based on indices
'''