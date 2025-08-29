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

ALGO_METADATA = {
    'momma': {
        'qvar':'Q'
    },
    'hivdi': {
        'qvar':'reach/Q'
    },
    'neobam':{
        'qvar':'q/q'
    },
    'metroman':{
        'qvar':'average/allq'
    },
    'sic4dvar':{
        'qvar':'Q_da'
    },
    'sad':{
        'qvar':'Qa'
    }
}

def process_reach(reach_id, mntdir):
  """Compute consensus for a single reach.

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
  for algo, metadata in ALGO_METADATA.items():
    infile = mntdir / 'flpe' / algo / f'{reach_id}_{algo}.nc'
    if not os.path.exists(infile):
      continue
    try:
      with Dataset(infile, 'r') as ds:
        arr = ds[metadata['qvar']][:].filled(np.nan)
        # treat negative discharge as NaN
        arr[arr<0] = np.nan
        # ignore algos with no nonnegative discharge
        if not np.any(arr>=0):
          continue
        arrs.append(arr)
        included_algos.append(algo)
    except (IOError, OSError):
      continue

  if not len(arrs):
    print(f"No data for reach '{reach_id}'")
    return

  consensus_arr = np.nanmedian(np.stack(arrs, axis=0), axis=0)
  
  outdir = mntdir / 'flpe' / 'consensus'
  if not os.path.exists(outdir):
    os.makedirs(outdir, exist_ok=True)

  outfile = outdir / f'{reach_id}_consensus.nc'

  with Dataset(outfile, 'w', format="NETCDF4") as dsout:
    dsout.n_algos = str(len(included_algos))
    dsout.contributing_algos = included_algos
    dsout.createDimension("nt", len(consensus_arr))
    consensus_q = dsout.createVariable("consensus_q", "f8", ("nt"), fill_value=np.nan)
    consensus_q.long_name = 'consensus discharge'
    consensus_q[:] = consensus_arr

def run_consensus(mntdir, indices):
  """Run consensus algorithm on a set of reaches.

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
