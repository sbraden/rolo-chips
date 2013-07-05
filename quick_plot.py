# python
# Sarah Braden

import pandas as pd
import numpy as np
import math
from matplotlib import rcParams 
from matplotlib import pyplot as plt
import itertools
import caltools # my module

# Gets rid of overlapping labels at the origin
rcParams['xtick.major.pad'] = 8
rcParams['ytick.major.pad'] = 8


def main():

    roloroot = '/home/sbraden/Dropbox/calibration/ROLO/datasetChips/'
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

    # Plotting
    index = 24 #695 nm 17 => 353 nm
    band = int(lookup.loc[index].values)
    plt.figure()
    chips[chips['filter_id'] == index].plot(x='topo_phase', y=feature, marker='o', lw=0)
    title = feature+'-phase-curve-ROLO-band-'+str(band)
    plt.title(title)
    plt.ylabel('I/F')
    plt.savefig(title+'.png', dpi=200)
    plt.show()
    plt.close

if __name__ == '__main__':
  main()