from __future__ import print_function
import sys, datetime
from multiprocessing import Pool
from clone_simulation_hpc2hpc import clone_simulation_hpc2hpc

from COMPS import Client
from COMPS.Data import Experiment, SimulationFile, QueryCriteria, Configuration, Simulation
from COMPS.Data.Simulation import SimulationState

import pdb

compshost = 'https://comps.idmod.org'


def should_rerun_sim(s):
    """
    Method to validate that this sim should be rerun (failed, is missing some
    expected output files, etc...
    """

    if s.state == SimulationState.Failed:
        return True

    try:
        Client.auth_manager()
    except:
        Client.login(compshost)

    fi = s.retrieve_output_file_info(None)

    if not any(filter(lambda x: x.path_from_root == "output" and x.friendly_name.startswith("MalariaSummaryReport_"), fi)):
        return True

    return False


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print('\r\nUsage:\r\n\t{0} (id-of-exp)'.format(sys.argv[0]))
        exit()

    Client.login(compshost)

    exp_id = sys.argv[1]
    exp = Experiment.get(exp_id)

    sims = exp.get_simulations()

    print(datetime.datetime.now())

    with Pool() as p:
        results = p.map(should_rerun_sim, sims)

    sims_to_rerun = [ sims[i] for i in filter(lambda x: results[x] == True, range(len(results))) ]

    print(datetime.datetime.now())

    if len(sims_to_rerun) == 0:
        print("Found no sims to rerun")
        exit(0)

    print("Found {} sims to rerun".format(len(sims_to_rerun)))

    expid = None

    for fs in sims_to_rerun:
        new_sim_id = clone_simulation_hpc2hpc(fs.id, expid)
        if expid == None:
            new_sim = Simulation.get(new_sim_id)
            expid = new_sim.experiment_id

    print("Recommissioning experiment")

    if exp.id == expid:
        exp.commission()

        # (my) existing exp... retag the sims I'm rerunning so I can delete them later
        for fs in sims_to_rerun:
            fs.merge_tags({'ClonedToRerun': None})
    else:
        exp2 = Experiment.get(expid)
        exp2.commission()

    print("Done")