import json
import os
import pandas as pd
import numpy as np

json_fname = "/Users/bertozzivill/Desktop/MalariaSummaryReport_10.json"
with open(json_fname) as f:
    report = json.loads(f.read())


timeinterval = report["Metadata"]["Reporting_Interval"]
agebins = report["Metadata"]["Age Bins"]
incdata = report["DataByTimeAndAgeBins"]["Annual Clinical Incidence by Age Bin"]
severedata =  report["DataByTimeAndAgeBins"]["Annual Severe Incidence by Age Bin"]
popdata =  report["DataByTimeAndAgeBins"]["Average Population by Age Bin"]
timedata = report["DataByTime"]["Time Of Report"]

max_years = int(max(timedata)/timeinterval)
incidence = [np.average(incdata[x], weights=popdata[x]) for x in range(0, max_years)]
severe_incidence = [severedata[x][0] for x in range(0, max_years)]


