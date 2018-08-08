from COMPS.Data import Simulation

from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser
from simtools.Utilities.COMPSUtilities import get_experiment_by_id, get_simulations_from_big_experiments

EXPERIMENT_ID = "3421c801-1b96-e811-a2c0-c4346bcb7275"
SIMULATIONS_PER_EXPERIMENT = 62500
SUITE_NAME = "my_suite"


if __name__ == "__main__":
    SetupParser.init("HPC")

    # Retrieve experiment
    e = get_experiment_by_id(EXPERIMENT_ID)
    print("The following experiment will be slitted:")
    print(e)

    # Retrieve simulations
    sims = get_simulations_from_big_experiments(e.id)
    print("{} simulations retrieved".format(len(sims)))

    # Calculate how many experiments we need
    exp_num = round(len(sims)/SIMULATIONS_PER_EXPERIMENT)
    print("Creating {} experiments with {} simulations each (max).".format(exp_num, SIMULATIONS_PER_EXPERIMENT))

    # Create the suite
    em = ExperimentManagerFactory.init()
    s = em.create_suite(SUITE_NAME)
    print("Suite named {} created with the id: {}".format(SUITE_NAME, s))

    # Create the experiments
    experiments = []
    for i in range(exp_num):
        em.create_experiment("{}_{}".format(SUITE_NAME, i), suite_id=s)
        experiments.append(em.comps_experiment)
        print("Created the experiment: {}".format(em.experiment.exp_name))

    exp_index = 0
    sims_in_exp = 0
    # For each simulations batch them
    for s in sims:
        s.experiment_id = experiments[exp_index].id
        sims_in_exp += 1

        if sims_in_exp == SIMULATIONS_PER_EXPERIMENT:
            print("Experiment {} full...".format(experiments[exp_index].name))
            Simulation.save_all(save_batch_callback=None, save_semaphore=Simulation.get_save_semaphore())
            exp_index += 1
            sims_in_exp = 0

    print("Processing done. Use the new suite id for analysis:")
    print(s)
