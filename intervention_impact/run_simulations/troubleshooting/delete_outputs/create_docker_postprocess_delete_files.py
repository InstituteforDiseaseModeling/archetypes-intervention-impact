import time, os, shutil, sys, uuid

from COMPS import Client
from COMPS.Data import WorkItem, WorkItemFile, QueryCriteria
from COMPS.Data.WorkItem import WorkItemState, WorkerOrPluginKey, RelationType

compshost = 'https://comps.idmod.org'
compsenv = 'Belegost'

if len(sys.argv) < 2 or len(sys.argv) > 3:
    print('\r\nUsage:\r\n\t{0} (comma-delimited-list-of-exp-ids-to-post-process) [0/1 to indicate overwrite existing postprocessing]'.format(sys.argv[0]))
    exit()

expids = [ uuid.UUID(eid) for eid in sys.argv[1].split(',') ]
overwrite = bool(int(sys.argv[2] if len(sys.argv) == 3 else 0))

Client.login(compshost)

# Create a work-item (locally)
wi = WorkItem('Postprocess DTK Simulations - delete files', WorkerOrPluginKey('DockerWorker', '1.0.0.0_RELEASE'), compsenv)

tags = {
        'WorkItem type': 'Docker',
        'overwrite': str(overwrite)
}

wi.set_tags(tags)

workorder_string = \
"""{
    "WorkItem_Type" : "DockerWorker",
    "Execution": {
        "ImageName": "ubuntu18.04_python3.6_dtk-tools1.1.6",
        "Command": "python3 postprocess_experiment_delete_files.py"
    }
}"""

wi.add_work_order(data=bytes(workorder_string, 'utf-8'))

additional_files = [
                       # filetype  filename   strip_carriage_returns
                       # ('input', 'run.sh',     1),
                        ('input', 'postprocess_experiment_delete_files.py')
                   ]

# add the linked files
for tup in additional_files:
    if len(tup) == 3 and tup[2] == 1:   # strip carriage returns
        with open(os.path.join('files', tup[1]), 'rb') as f:
            wi.add_file(WorkItemFile(tup[1], tup[0], ''), data=f.read().replace('\r\n', '\n'))
    else:
        wi.add_file(WorkItemFile(tup[1], tup[0], ''), file_path=tup[1])

# Save the work-item to the server
wi.save()

for eid in expids:
    print('Adding relationship to experiment {0}'.format(str(eid)))
    wi.add_related_experiment(eid, RelationType.DependsOn)

wi.refresh()

print('Created work-item {0}'.format(wi.id))
print('Commissioning...')

wi.commission()

print('Commissioned')
print('Refreshing WorkItem state until it completes')

print('State -> {}'.format(wi.state.name))

while wi.state not in [WorkItemState.Succeeded, WorkItemState.Failed, WorkItemState.Canceled]:
    time.sleep(5)
    wi.refresh()
    print('State -> {}'.format(wi.state.name))

print('Worker complete: ' + wi.state.name)

if wi.state is not WorkItemState.Succeeded:
    exit(-1)

print('Showing stdout of workitem:')

barr_out = wi.retrieve_output_files(['stdout.txt'])

for line in barr_out[0].decode('utf-8').splitlines():
    print(' > {0}'.format(line))