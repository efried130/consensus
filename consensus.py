"""Module to compute consensus from six base FLPEs

Runs on a single reach or set of reaches and requires JSON data for reach retrieved by
AWS Batch index.
"""

import argparse as ap
from pathlib import Path
import json
from netCDF4 import Dataset
import numpy as np
import os
import datetime
import pandas as pd

ALGO_METADATA = {
    'momma': {
        'qvar':'Q',
        'time':'time_str'
    },
    'hivdi': {
        'qvar':'reach/Q',
        'time':'time'
    },
    'neobam':{
        'qvar':'q/q',
        'time':'time_str'
    },
    'metroman':{
        'qvar':'average/allq',
        'time':'time_str'
    },
    'sic4dvar':{
        'qvar':'Q_da',
        'time':'times'
    },
    'sad':{
        'qvar':'Qa',
        'time':'time_str'
    }
}

FILL_VALUE = -999999999999.0
FILL_VALUE_STR = "no_data"

def remove_low_cv_and_recalc_consensus(arrs, time_arrs, CV_thresh, included_algos):
    """
    For a list of discharge arrays:
    - Removes arrays with CV < threshold
    - Recalculates consensus using the remaining arrays
    Parameters
    ----------
    arrs : list of np.ndarray
        Discharge arrays from each algorithm.
    CV_thresh : float
        Coefficient of variation threshold below which arrays are excluded.

    Returns
    -------
    np.ndarray
        Cleaned and recalculated consensus array.
    """
    
    cv_arrs = []
    cv_included_algos = []
    cv_time_arrs = []
    
    for i, arr in enumerate(arrs):
        mean = np.nanmean(arr)
        std = np.nanstd(arr)
        cv = std / mean if mean != 0 else np.nan
        
        if not np.isnan(cv) and cv > CV_thresh:
            cv_arrs.append(arr)
            cv_included_algos.append(included_algos[i])
            cv_time_arrs.append(time_arrs[i])

    if not len(cv_arrs):
        print("All algorithms removed due to low CV; returning NaN array and no included algos.")
        return np.full_like(arrs[0], np.nan), np.full_like(arrs[0], "no_data", dtype=object), []

    # Compute median consensus
    consensus_arr = np.nanmedian(np.stack(cv_arrs, axis=0), axis=0)
    selected_time_arr = time_arrs[0]


    # # Debug print
    # print(f"Selected algo(s) for consensus: {cv_included_algos}")
    # print(f"Time array length: {len(selected_time_arr)}, Consensus array shape: {cv_arrs[0].shape}")
    # print(f"Time array : {(selected_time_arr)}, Consensus array : {cv_arrs}")

    return consensus_arr, selected_time_arr, cv_included_algos

def process_reach(reach_id, mntdir):
    """
    Compute consensus for a single reach.

    Parameters
    ----------
    mntdir: Path
    path to base mount directory
    reach_id: int
    ID of reach to process
    """

    print('reach', reach_id)
    included_algos = []
    arrs = []
    time_arrs = []
    for algo, metadata in ALGO_METADATA.items():
        infile = mntdir / 'flpe' / algo / f'{reach_id}_{algo}.nc'
        if not os.path.exists(infile):
            continue
        try:
            with Dataset(infile, 'r') as ds:
                arr = ds[metadata['qvar']][:].filled(np.nan)
                
                algo_time = ds.variables[metadata['time']][:]
                if algo == 'sic4dvar':
                    
                    mask = np.ma.getmaskarray(algo_time)
                    valid_indexes = [i for i in range(algo_time.shape[0])]

                    if valid_indexes:
                        valid_sic_str = [algo_time[i] for i in valid_indexes]  # keep all for indexing
                        swot_ts = datetime.datetime(2000,1,1,0,0,0)

                        # convert masked values to np.nan
                        valid_sic_str = np.array(valid_sic_str, dtype=float)
                        valid_sic_str = np.where(np.ma.getmaskarray(valid_sic_str), np.nan, valid_sic_str)

                        algo_time = np.array([
                            (swot_ts + datetime.timedelta(days=t)).strftime("%Y-%m-%dT%H:%M:%SZ") if not np.isnan(t) else None
                            for t in valid_sic_str
                        ])
                
                time=algo_time
                
                # treat negative discharge as NaN
                arr[arr<0] = np.nan
                # ignore algos with no nonnegative discharge
                if not np.any(arr>=0):
                    continue
                    return consensus_arr, consensus_time_arr, []


                arrs.append(arr)
                time_arrs.append(time)
                included_algos.append(algo)
                 
        except (IOError, OSError):
            continue               
    if not len(arrs):
        print(f"No data for reach '{reach_id}'")
        return


    ##CHOOSE WHETHER TO APPLY CV FILTER HERE
    # consensus_arr, time_arr = np.nanmedian(np.stack(arrs, axis=0), axis=0), time_arrs[0]
    consensus_arr, time_arr, included_algos = remove_low_cv_and_recalc_consensus(arrs=arrs, time_arrs=time_arrs, CV_thresh=0.5, included_algos=included_algos)

    
    #Build nc file
    outdir = mntdir / 'flpe' / 'consensus'
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    outfile = outdir / f'{reach_id}_consensus.nc'

    with Dataset(outfile, 'w', format="NETCDF4") as dsout:
        dsout.n_algos = str(len(included_algos))
        dsout.contributing_algos = included_algos
        
        #Add consensus Q
        dsout.createDimension("nt", len(consensus_arr))
        consensus_q = dsout.createVariable("consensus_q", "f8", ("nt"), fill_value=FILL_VALUE)
        consensus_q.long_name = 'consensus discharge'
        consensus_q.short_name = "discharge_consensus"
        consensus_q.tag_basic_expert = "Basic"
        consensus_q.units = "m^3/s"
        #consensus_q.quality_flag = "dschg_c_q" Doesn't yet exist
        consensus_q.valid_min = -10000000.0
        consensus_q.valid_max = 10000000.0
        # consensus_q.coordinates = "p_lon p_lat"
        consensus_q.comment = "Discharge from the consensus discharge algorithm."
        
        #Add consensus time_str
        consensus_time_str = dsout.createVariable("time_str", str, ("nt"), fill_value="no_data")
        consensus_time_str.long_name = "time (UTC)"
        consensus_time_str.standard_name = "time"
        consensus_time_str.short_name = "time_string"
        consensus_time_str.calendar = "gregorian"
        consensus_time_str.tag_basic_expert = "Basic"
        consensus_time_str.comment = (
            "Time string giving UTC time. The format is YYYY-MM-DDThh:mm:ssZ, "
            "where the Z suffix indicates UTC time."
        )

        #Fill as needed
        consensus_arr_filled = np.where(np.isnan(consensus_arr), FILL_VALUE, consensus_arr)
        time_arr_filled = [t if t is not None else FILL_VALUE_STR for t in time_arr]

        # Write values
        consensus_q[:] = consensus_arr_filled
        consensus_time_str[:] = np.array(time_arr_filled, dtype="O")
        
def run_consensus(mntdir, indices):
    """
    Run consensus algorithm on a set of reaches.

    Parameters
    ----------
    mntdir: Path
    path to base mount directory
    indices: list
    offsets of reaches to process
    """

    reachfile = mntdir / 'input' / 'reaches.json'

    with open(reachfile, 'r') as fp:
        reaches = json.load(fp)
        reach_ids = [reaches[i]['reach_id'] for i in indices]

    for reach_id in reach_ids:
        process_reach(reach_id, mntdir)

def parse_range(index_str):
    """Parse a range string into a list of integers."""

    indices = []
    try:
        for part in index_str.strip().split(","):
            part = part.strip()
            if "-" in part:
                start, end = part.split("-")
                indices.extend(list(range(int(start), int(end)+1)))
            else:
                indices.append(int(part))
    except (IndexError, ValueError, TypeError):
        print(f"cannot parse range string: '{index_str}'. Must be either a single integer, a range such as 1-100, or a comma separated lists of integers and/or ranges")
    
    # remove duplicates and sort
    return sorted(list(set(indices)))

if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument("--mntdir", type=str, default="/mnt", help="Mount directory.")
    parser.add_argument("-i", "--index", type=parse_range, required=True)
    args = parser.parse_args()
    run_consensus(Path(args.mntdir), args.index)
