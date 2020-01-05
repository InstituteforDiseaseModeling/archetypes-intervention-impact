from __future__ import print_function
import sys

from COMPS import Client
from COMPS.Data import Experiment, SimulationFile, QueryCriteria, Configuration, Simulation

compshost = 'https://comps.idmod.org'
state_to_delete = 'Failed'

if len(sys.argv) < 2 or len(sys.argv) > 3:
    print('\r\nUsage:\r\n\t{0} (id-of-exp) [state-to-delete]'.format(sys.argv[0]))
    exit()

Client.login(compshost)

exp = Experiment.get(sys.argv[1])
if len(sys.argv) == 3:
    state_to_delete = sys.argv[2]

sims_to_delete = exp.get_simulations(QueryCriteria().select('id').where('state='+state_to_delete))

if len(sims_to_delete) == 0:
    print("Found no {0} sims to delete".format(state_to_delete))
    exit(0)

for fs in sims_to_delete:
    fs.delete()

print("Done")