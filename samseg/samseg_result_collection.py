import os
import dill
import pandas as pd

# the purpose of this script is to go through all the directories created for/by SAMSEG and to assemble the segmentation results
# this will have to be done in a different manner for subjects with different sets of brainscans (with/without follow-up, with/without FLAIR, etc.)

# load the data prepared with slurm_samseg_creator_preparation.py
dill.load_session("/data/pt_life_freesurfer/samseg/scripts/prepared_vars.pkl")

# let's start with subjects that only have one timepoint and both a T1 and a FLAIR
# load a subject to get the names of all measured regions to use as column names for the summary table
df = pd.read_csv('/data/pt_life_freesurfer/samseg/DE12F2B3CF_bl/samseg.stats', names = ["region", "size", "unit"])
# remove "# Measure " from the region labels
df['region'] = df['region'].str[10:]
# get the names of the segmented brain regios to use as cols
brain_regions = df['region'].tolist()

# create a dataframe to collect all segmentation results
samseg_results = pd.DataFrame({
    'participant': [None] * 3729, # 3729 is the total number of timepoints over all particpants
    'timepoint': [None] * 3729,
    'longitudinal': [None] * 3729,
    'T1': [None] * 3729,
    'FLAIR': [None] * 3729,
    'sbTIV': [None] * 3729
})

for col in brain_regions:
    samseg_results[col] = [None] * 3729
# fill the prepared dataframe with the basic information on this group of participants
length_onetpboth = len(onetpboth)
samseg_results.iloc[:length_onetpboth, 2] = 0  # assign the value 0 to the "longitudinal" column, only one timepoint, no longitudinal pipeline
samseg_results.iloc[:length_onetpboth, 3] = 1  # assign the value 1 to the "T1" column, T1 was available
samseg_results.iloc[:length_onetpboth, 4] = 1  # assign the value 1 to the "FLAIR" column, FLAIR was available
# go through all participants in the list onetpboth and extract the samseg results
# initiate counter to assign proper row
n = 0
for participant in onetpboth:
    samseg_results.iloc[n,0] = participant[4:] # assign participant's name to the row of the samseg_results dataframe
    # check if the single available timepoint is a baseline or follow-up measurement
    if os.path.exists("/data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/"):
        samseg_results.iloc[n, 1] = "BL" # assign the value "BL" to the timepoint col
        sbtiv = pd.read_csv(filepath_or_buffer= "/data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/sbtiv.stats", names = ["region", "size", "unit"])
        samseg_results.iloc[n, 5] = sbtiv.iloc[0,1] # get the ICV from the sbtiv.stats file
        indiv_samseg = pd.read_csv(filepath_or_buffer= "/data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/samseg.stats", names = ["region", "size", "unit"]) # read in the particpant's samseg results
        samseg_results.iloc[n,6:] = indiv_samseg["size"] # assign the participant's regions' sizes to the columns in the participant's row in the summary table
    else:
        samseg_results.iloc[n, 1] = "FU"  # assign the value "FU" to the timepoint col
        sbtiv = pd.read_csv(filepath_or_buffer="/data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/sbtiv.stats",
                            names=["region", "size", "unit"])
        samseg_results.iloc[n, 5] = sbtiv.iloc[0,1]  # get the ICV from the sbtiv.stats file
        indiv_samseg = pd.read_csv(
            filepath_or_buffer="/data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/samseg.stats",
            names=["region", "size", "unit"])  # read in the particpant's samseg results
        samseg_results.iloc[n, 6:] = indiv_samseg["size"]  # assign the participant's regions' sizes to the columns in the participant's row in the summary table
    n += 1 # update the counter
# now let's continue with the subjects that have one timepoint and only a T1 and no FLAIR scans
length_onetpt1 = len(onetpt1)
samseg_results.iloc[length_onetpboth:length_onetpboth+length_onetpt1, 2] = 0  # assign the value 0 to the "longitudinal" column, only one timepoint, no longitudinal pipeline
samseg_results.iloc[length_onetpboth:length_onetpboth+length_onetpt1, 3] = 1  # assign the value 1 to the "T1" column, T1 was available
samseg_results.iloc[length_onetpboth:length_onetpboth+length_onetpt1, 4] = 0  # assign the value 0 to the "FLAIR" column, FLAIR wasn't available
# go through all participants in the list onetpt1 and extract the samseg results
# continue to use the counter to assign the proper row
for participant in onetpt1:
    samseg_results.iloc[n,0] = participant[4:] # assign participant's name the row of the samseg_results dataframe
    # check if the single available timepoint is a baseline or follow-up measurement
    if os.path.exists("/data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/"):
        samseg_results.iloc[n, 1] = "BL" # assign the value "BL" to the timepoint col
        sbtiv = pd.read_csv(filepath_or_buffer="/data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/sbtiv.stats",
                            names=["region", "size", "unit"])
        samseg_results.iloc[n, 5] = sbtiv.iloc[0,1]  # get the ICV from the sbtiv.stats file
        indiv_samseg = pd.read_csv(filepath_or_buffer= "/data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/samseg.stats", names = ["region", "size", "unit"]) # read in the particpant's samseg results
        # participants without FLAIR have different cols in the samseg.stats file
        # they do not have a lesions or a left undetermined col (our last two cols)
        # instead they have an additional WM-hypointensities (33rd) and non-WM-hypointensities (43rd) col which we will ignore
        relevant_cols = list(indiv_samseg["size"][0:32]) + list(indiv_samseg["size"][33:42])
        relevant_cols.append(float(indiv_samseg["size"][43]))
        samseg_results.iloc[n,6:48] = relevant_cols # assign the participant's regions' sizes to the columns in the participant's row in the summary table, leave last 2 cols empty
    else:
        samseg_results.iloc[n, 1] = "FU"  # assign the value "FU" to the timepoint col
        sbtiv = pd.read_csv(filepath_or_buffer="/data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/sbtiv.stats",
                            names=["region", "size", "unit"])
        samseg_results.iloc[n, 5] = sbtiv.iloc[0,1]  # get the ICV from the sbtiv.stats file
        indiv_samseg = pd.read_csv(
            filepath_or_buffer="/data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/samseg.stats",
            names=["region", "size", "unit"])  # read in the particpant's samseg results
        # participants without FLAIR have different cols in the samseg.stats file
        # they do not have a lesions or a left undetermined col (our last two cols)
        # instead they have an additional WM-hypointensities (33rd) and non-WM-hypointensities (43rd) col which we will ignore
        relevant_cols = list(indiv_samseg["size"][0:32]) + list(indiv_samseg["size"][33:42])
        relevant_cols.append(float(indiv_samseg["size"][43]))
        samseg_results.iloc[n,6:48] = relevant_cols  # assign the participant's regions' sizes to the columns in the participant's row in the summary table, leave last 2 cols empty
    n += 1 # update the counter

# let's now continue with subjects with two available timepoints and both modalities available at both timepoints
length_twotpbothboth = len(twotpbothboth)*2 # double it as there are 2 tps per participant
samseg_results.iloc[length_onetpboth+length_onetpt1:length_onetpboth+length_onetpt1+length_twotpbothboth, 2] = 1  # assign the value 1 to the "longitudinal" column, the longitudinal pipeline was used
samseg_results.iloc[length_onetpboth+length_onetpt1:length_onetpboth+length_onetpt1+length_twotpbothboth, 3] = 1  # assign the value 1 to the "T1" column, T1 was available
samseg_results.iloc[length_onetpboth+length_onetpt1:length_onetpboth+length_onetpt1+length_twotpbothboth, 4] = 1  # assign the value 1 to the "FLAIR" column, FLAIR was available
# go through all participants in the list twotpbothboth and extract the samseg results
# continue to use the counter to assign the proper row
for participant in twotpbothboth:
    samseg_results.iloc[n,0] = participant[4:] # assign participant's name the row of the samseg_results dataframe for the BL measurement
    samseg_results.iloc[n, 1] = "BL" # assign the value "BL" to the timepoint col
    if not os.path.exists("/data/pt_life_freesurfer/samseg/" + participant[4:] + "/tp001/samseg.stats"):
        if os.path.exists("/data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/sbtiv.stats"):
            sbtiv = pd.read_csv(filepath_or_buffer="/data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/sbtiv.stats",
                names=["region", "size", "unit"])
            samseg_results.iloc[n, 5] = sbtiv.iloc[0,1]  # get the ICV from the sbtiv.stats file
            indiv_samseg = pd.read_csv(filepath_or_buffer="/data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/samseg.stats",
                names=["region", "size", "unit"])  # read in the particpant's samseg results
            samseg_results.iloc[n, 6:] = indiv_samseg["size"]  # assign the participant's regions' sizes to the columns in the participant's row in the summary table
        else:
        samseg_results.iloc[n,5:] = None
    else:
        sbtiv = pd.read_csv(filepath_or_buffer="/data/pt_life_freesurfer/samseg/" + participant[4:] + "/tp001/sbtiv.stats",
                            names=["region", "size", "unit"])
        samseg_results.iloc[n, 5] = sbtiv.iloc[0,1]  # get the ICV from the sbtiv.stats file
        indiv_samseg = pd.read_csv(filepath_or_buffer= "/data/pt_life_freesurfer/samseg/" + participant[4:] + "/tp001/samseg.stats", names = ["region", "size", "unit"]) # read in the particpant's samseg results
        samseg_results.iloc[n,6:] = indiv_samseg["size"] # assign the participant's regions' sizes to the columns in the participant's row in the summary table
    n += 1 # update the counter
    samseg_results.iloc[n, 0] = participant[4:]
    samseg_results.iloc[n, 1] = "FU"  # assign the value "FU" to the timepoint col
    if not os.path.exists("/data/pt_life_freesurfer/samseg/" + participant[4:] + "/tp002/samseg.stats"):
        if os.path.exists("/data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/sbtiv.stats"): # check if cross-sectional results exist if the longitudinal pipeline failed
            sbtiv = pd.read_csv(
                filepath_or_buffer="/data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/sbtiv.stats",
                names=["region", "size", "unit"])
            samseg_results.iloc[n, 5] = sbtiv.iloc[0, 1]  # get the ICV from the sbtiv.stats file
            indiv_samseg = pd.read_csv(
                filepath_or_buffer="/data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/samseg.stats",
                names=["region", "size", "unit"])  # read in the particpant's samseg results
            samseg_results.iloc[n, 6:] = indiv_samseg[
                "size"]  # assign the participant's regions' sizes to the columns in the participant's row in the summary table
        else:
            samseg_results.iloc[n, 5:] = None
    else:
        sbtiv = pd.read_csv(
            filepath_or_buffer="/data/pt_life_freesurfer/samseg/" + participant[4:] + "/tp002/sbtiv.stats",
            names=["region", "size", "unit"])
        samseg_results.iloc[n, 5] = sbtiv.iloc[0,1]  # get the ICV from the sbtiv.stats file
        indiv_samseg = pd.read_csv(
        filepath_or_buffer="/data/pt_life_freesurfer/samseg/" + participant[4:] + "/tp002/samseg.stats",
        names=["region", "size", "unit"])  # read in the particpant's samseg results
        samseg_results.iloc[n, 6:] = indiv_samseg["size"]  # assign the participant's regions' sizes to the columns in the participant's row in the summary table
    n += 1 # update the counter

# let's continue with the subject that had both modalities at baseline and only T1 and FU
length_twotpbotht1 = len(twotpbotht1)*2 # double it as there are 2 tps per participant
samseg_results.iloc[length_onetpboth+length_onetpt1+length_twotpbothboth:length_onetpboth+length_onetpt1+length_twotpbothboth+length_twotpbotht1, 2] = 0  # assign the value 0 to the "longitudinal" column, the longitudinal pipeline was not used due to differing modalities at the two timepoints
samseg_results.iloc[length_onetpboth+length_onetpt1+length_twotpbothboth:length_onetpboth+length_onetpt1+length_twotpbothboth+length_twotpbotht1, 3] = 1  # assign the value 1 to the "T1" column, T1 was available
samseg_results.iloc[length_onetpboth+length_onetpt1+length_twotpbothboth:length_onetpboth+length_onetpt1+length_twotpbothboth+1, 4] = 1  # assign the value 1 to the "FLAIR" column, FLAIR was available for the baseline
samseg_results.iloc[length_onetpboth+length_onetpt1+length_twotpbothboth+1:length_onetpboth+length_onetpt1+length_twotpbothboth+2, 4] = 0  # assign the value 0 to the "FLAIR" column, FLAIR wasn't available for the follow-up

for participant in twotpbotht1:
    samseg_results.iloc[n,0] = participant[4:] # assign participant's name the row of the samseg_results dataframe for the BL measurement
    samseg_results.iloc[n, 1] = "BL" # assign the value "BL" to the timepoint col
    sbtiv = pd.read_csv(filepath_or_buffer="/data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/sbtiv.stats",
                        names=["region", "size", "unit"])
    samseg_results.iloc[n, 5] = sbtiv.iloc[0,1]  # get the ICV from the sbtiv.stats file
    indiv_samseg = pd.read_csv(filepath_or_buffer= "/data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/samseg.stats", names = ["region", "size", "unit"]) # read in the particpant's samseg results
    samseg_results.iloc[n,6:] = indiv_samseg["size"] # assign the participant's regions' sizes to the columns in the participant's row in the summary table
    n += 1 # update the counter
    samseg_results.iloc[n, 0] = participant[4:]
    samseg_results.iloc[n, 1] = "FU"  # assign the value "FU" to the timepoint col
    sbtiv = pd.read_csv(filepath_or_buffer="/data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/sbtiv.stats",
                        names=["region", "size", "unit"])
    samseg_results.iloc[n, 5] = sbtiv.iloc[0,1]  # get the ICV from the sbtiv.stats file
    indiv_samseg = pd.read_csv(
        filepath_or_buffer="/data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/samseg.stats",
        names=["region", "size", "unit"])  # read in the particpant's samseg results
    # participants without FLAIR have different cols in the samseg.stats file
    # they do not have a lesions or a left undetermined col (our last two cols)
    # instead they have an additional WM-hypointensities (33rd) and non-WM-hypointensities (43rd) col which we will ignore
    relevant_cols = list(indiv_samseg["size"][0:32]) + list(indiv_samseg["size"][33:42])
    relevant_cols.append(float(indiv_samseg["size"][43]))
    samseg_results.iloc[n,6:48] = relevant_cols  # assign the participant's regions' sizes to the columns in the participant's row in the summary table, leave last 2 cols empty
    n += 1 # update the counter


# let's finish it off with the subject that had both modalities at FU and only T1 at baseline
length_twotpt1both = len(twotpt1both)*2 # double it as there are 2 tps per participant
samseg_results.iloc[length_onetpboth+length_onetpt1+length_twotpbothboth+length_twotpbotht1:length_onetpboth+length_onetpt1+length_twotpbothboth+length_twotpbotht1+length_twotpt1both, 2] = 0  # assign the value 0 to the "longitudinal" column, the longitudinal pipeline was not used due to differing modalities at the two timepoints
samseg_results.iloc[length_onetpboth+length_onetpt1+length_twotpbothboth+length_twotpbotht1:length_onetpboth+length_onetpt1+length_twotpbothboth+length_twotpbotht1+length_twotpt1both, 3] = 1  # assign the value 1 to the "T1" column, T1 was available
samseg_results.iloc[length_onetpboth+length_onetpt1+length_twotpbothboth+length_twotpbotht1:length_onetpboth+length_onetpt1+length_twotpbothboth+length_twotpbotht1+1, 4] = 0  # assign the value 0 to the "FLAIR" column, FLAIR was not available for the baseline
samseg_results.iloc[length_onetpboth+length_onetpt1+length_twotpbothboth+length_twotpbotht1+1:length_onetpboth+length_onetpt1+length_twotpbothboth+length_twotpbotht1+2, 4] = 1  # assign the value 0 to the "FLAIR" column, FLAIR was available for the follow-up

for participant in twotpt1both:
    samseg_results.iloc[n,0] = participant[4:] # assign participant's name the row of the samseg_results dataframe for the BL measurement
    samseg_results.iloc[n, 1] = "BL" # assign the value "BL" to the timepoint col
    sbtiv = pd.read_csv(filepath_or_buffer="/data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/sbtiv.stats",
                        names=["region", "size", "unit"])
    samseg_results.iloc[n, 5] = sbtiv.iloc[0,1]  # get the ICV from the sbtiv.stats file
    indiv_samseg = pd.read_csv(filepath_or_buffer= "/data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/samseg.stats", names = ["region", "size", "unit"]) # read in the particpant's samseg results
    # participants without FLAIR have different cols in the samseg.stats file
    # they do not have a lesions or a left undetermined col (our last two cols)
    # instead they have an additional WM-hypointensities (33rd) and non-WM-hypointensities (43rd) col which we will ignore
    relevant_cols = list(indiv_samseg["size"][0:32]) + list(indiv_samseg["size"][33:42])
    relevant_cols.append(float(indiv_samseg["size"][43]))
    samseg_results.iloc[n,6:48] = relevant_cols  # assign the participant's regions' sizes to the columns in the participant's row in the summary table, leave last 2 cols empty
    n += 1 # update the counter
    samseg_results.iloc[n, 0] = participant[4:]
    samseg_results.iloc[n, 1] = "FU"  # assign the value "FU" to the timepoint col
    sbtiv = pd.read_csv(filepath_or_buffer="/data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/sbtiv.stats",
                        names=["region", "size", "unit"])
    samseg_results.iloc[n, 5] = sbtiv.iloc[0,1]  # get the ICV from the sbtiv.stats file
    indiv_samseg = pd.read_csv(
        filepath_or_buffer="/data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/samseg.stats",
        names=["region", "size", "unit"])  # read in the particpant's samseg results
    samseg_results.iloc[n, 6:] = indiv_samseg["size"]  # assign the participant's regions' sizes to the columns in the participant's row in the summary table
    n += 1 # update the counter


samseg_results.to_csv('/data/pt_life_freesurfer/samseg/results_summary.csv')