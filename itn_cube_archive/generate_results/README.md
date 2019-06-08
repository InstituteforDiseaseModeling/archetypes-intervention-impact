Code for generating ITN coverage "cubes" originally by Sam Bhatt, refactored by Amelia Bertozzi-Villa. 

* sam_orig is the code directly from sam, unchanged except for some comments. 
* run_on_gcloud is a minimal refactoring of sam's code to run on the google cloud framework-- mostly filepath updates and some function refactoring to get the process running
* amelia_refactor is a complete refactor of the process, but preserving all bugs such that it replicates the outputs of run_on_glcoud exactly. 

The current amelia_refactor replicates the run_on_glcoud results exactly, including some suspected bugs, and preserves the original dataset formatting. I'm freezing this part of the repo now and moving the amelia_refactor folder to the map-itn-cube for bug fixes and further updates.

Last Updated: ABV, 8 June 2019
