import uuid
import time
import sys
import os
import re
import subprocess
from multiprocessing import Pool

from COMPS import Client
from COMPS.Data import Simulation, Configuration, Priority, QueryCriteria, SimulationFile, Experiment, WorkItem
from COMPS.Data.Simulation import SimulationState

compshost = 'https://comps.idmod.org'

files_to_delete = [
    'StdOut.txt'
]

def should_postprocess_sim(s, overwrite_existing=False):
    """
    Method to validate that this sim should be post-processed (succeeded/failed, had the correct
    files, was/wasn't already preprocessed, etc...
    """

    if s.state == SimulationState.Succeeded:
        return True

    return False

def postprocess_sim(s):
    """
    Controls what actual commands are run to preprocess the existing sim
    """

    s.refresh(QueryCriteria().select_children('hpc_jobs'))

    simpath = unc_path_to_docker_path(s.hpc_jobs[-1].working_directory)

    num_files_deleted = 0

    for f in files_to_delete:
        fp = os.path.join(simpath, f)
        if os.path.exists(fp):
            try:
                os.remove(fp)
                num_files_deleted += 1
            except:
                pass

    return (s, num_files_deleted)


#
# Probably no need to change anything below here...
#


def unc_path_to_docker_path(p):
    mapping = os.environ['COMPS_DATA_MAPPING'].split(';')
    return re.sub(mapping[1].replace('\\', '\\\\'), mapping[0], p, flags=re.IGNORECASE).replace('\\', '/')


if __name__ == "__main__":
    if len(sys.argv) > 2:
        print('\r\nUsage:\r\n\t{0} [comma-delimited-list-of-exp-ids-to-post-process]'.format(sys.argv[0]))
        exit()

    Client.login(compshost)

    wi = WorkItem.get(os.environ['COMPS_WORKITEM_GUID'], QueryCriteria().select_children(['tags']))

    if len(sys.argv) == 2:
        exps = [ Experiment.get(eid) for eid in sys.argv[1].split(',') ]
    else:
        exps = wi.get_related_experiments()

    for exp in exps:
        print('Experiment {0}'.format(str(exp.id)))

        sims = exp.get_simulations()

        filtered_sims = [ s for s in sims if should_postprocess_sim(s, wi.tags.get("overwrite") == "True") ]

        print("Found {} sims to post-process".format(len(filtered_sims)))

        with Pool() as p:
            results = p.map(postprocess_sim, filtered_sims)

        for s, result in results:
            print()
            print("Post-processed simulation {}, deleted {} files".format(str(s.id), result))
