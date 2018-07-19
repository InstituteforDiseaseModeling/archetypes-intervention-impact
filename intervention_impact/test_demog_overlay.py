import pandas as pd


from input_file_generation.add_properties_to_demographics import generate_demographics_properties, check_df_valid

# base demographics file
refdemo_fname = "sites/aba/demographics_aba.json"

as_overlay = True # generate demographics overlay or directly integrate into base demographics file
output_filename = "sites/aba/demographics_overlay_aba.json"

# read properties to overlay from csv or enter in dict (prop_by_node_dict)
# csv must be in same format as this example dictionary
# if not setting properties by node, pass empty dict.
# expected format for df:
# columns:
# 'node' : node id
# 'Property' : name of property to set
# 'Property_Type' : 'NP' or 'IP'
# 'Property_Value' : value of property to be set
# 'Initial_Distribution' : initial fraction of population to have the row's IP value. Ignored for NPs.
# for example
# prop_by_node_dict = { 'node' : [3, 3, 3, 3, 3],
#                       'Property' : ['NodeType', 'ForestGoing', 'ForestGoing', 'DrugStatus', 'DrugStatus'],
#                       'Property_Type' : ['NP', 'IP', 'IP', 'IP', 'IP'],
#                       'Property_Value' : ['Forest', 'LovesForest', 'HatesForest', 'None', 'RecentDrug'],
#                       'Initial_Distribution' : [1, 0.9, 0.1, 0.3, 0.7]}
from_csv = False
csv_fname = ''
prop_by_node_dict = { 'node' : [1475161276, 1475161276],
                      'Property' : ['NetUsage', 'NetUsage'],
                      'Property_Type' : ['IP', 'IP'],
                      'Property_Value' : [ 'HatesNets',
                                           'LovesNets'],
                      'Initial_Distribution' : [0.2, 0.8]}

# if NPs and IPs do not already exist in base demographics, specify them here; otherwise set to empty lists
IPs = [
    { 'Property' : 'NetUsage',
      'Values' : [ 'HatesNets',
                   'LovesNets'],
      'Initial_Distribution' : [0, 1],
      'Transitions' : [] }
]
NPs = [
]

df = pd.DataFrame(prop_by_node_dict)

if (not df.empty) and (not check_df_valid(df, IPs, NPs)):
    print('properties by node df is invalid')
    exit()
generate_demographics_properties(refdemo_fname, output_filename, as_overlay, IPs, NPs, df)

