# this is a script to prepare information for the creation of slurm scipts
# slurm_samseg_creator.py will create multiple slurm scripts to run freesurfer's samseg across all life subjects
# to be run parallelized via the slurm compute cluster
# if participants have both a T1 scan and a FLAIR scan both will be used
# if participants underwent MRIs at baseline and follow-up, the longitudinal version of samseg will be run
# in the longitudinal version a template based on all scans has to be created before the individual timepoints can be analysed
import os
from collections import OrderedDict
import dill

# we will create many slurm scripts
# the following function will create a header for the scripts and help keep this script leaner
# n - how many participants should be processed via the script
# m - starting with which participant from the current list
# c - how many cores should be requested
# name - name of the script, e.g. "slurm_onetp_coregistration_"
def slurmheader(c, n, m, name, mem = "12"):
      print("#!/bin/bash",
            "#SBATCH --job-name=" + name + str(m) + "_to_" + str(m+n-1),         # create a short name for your job",
            "#SBATCH -c " + str(c) + "                # core count",
            "#SBATCH --ntasks=" + str(c) + "              # total number of tasks across all nodes",
            "#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)",
            "#SBATCH --mem-per-cpu=" + mem + "G         # 12G to avoid error: processes may have been killed by the cgroup out-of-memory handler",
            "#SBATCH --time=100:00:00          # total run time limit (HH:MM:SS)",
            "#SBATCH --mail-type=begin        # send email when job begins",
            "#SBATCH --mail-type=end          # send email when job ends",
            "#SBATCH --mail-user=lammer@cbs.mpg.de",
            "",
            sep="\n",
            file=open(("/data/pt_life_freesurfer/samseg/scripts/" + name + str(m) + "_to_" + str(
                  m + n - 1)), "a"))

#specify the parent directory
parent_directory = "/data/p_life_raw/bids/"

#initialize the dictionary
participant_dict = {}

#walk through the directory structure
for folder in os.listdir(parent_directory):
    folder_path = os.path.join(parent_directory, folder)
    if os.path.isdir(folder_path):  #check if it's a directory
        subfolders = [subfolder for subfolder in os.listdir(folder_path) if os.path.isdir(os.path.join(folder_path, subfolder))]
        participant_dict[folder] = subfolders

# only use the 2nd measurement if participants were scanned twice at the same timepoint
for participant, sessions in participant_dict.items():
      if "ses-bl" and "ses-bl2" in sessions:
            sessions.remove("ses-bl")
      if "ses-fu" and "ses-fu2" in sessions:
            sessions.remove("ses-fu")

# turn it into an ordered directory with NA
participants = OrderedDict()
for participant, sessions in participant_dict.items():
      if "ses-bl" in sessions:
            participants[participant] = ["ses-bl"]
      elif "ses-bl2" in sessions:
            participants[participant] = ["ses-bl2"]
      else:
            participants[participant] = ["NA"]
      if "ses-fu" in sessions:
            participants[participant].append("ses-fu")
      elif "ses-fu2" in sessions:
            participants[participant].append("ses-fu2")
      else:
            participants[participant].append("NA")

# there are 2666 participants
# all have at least one timepoint
# 7 have NA for BL
# 1589 have NA for FU - ergo 1077 have FU

# to keep things tidy we will define paths as variables
# the following function does that for a participant
def pathtovariable(participant, sessions):
      dir_bl = "/data/p_life_raw/bids/" + participant + "/" + sessions[0] + "/"
      dir_fu = "/data/p_life_raw/bids/" + participant + "/" + sessions[1] + "/"
      # cave: there are some inconsistencies in the naming of t1 and FLAIR files at baseline
      t1_bl = ""
      for suffix in ["_ses-bl_acq-ADNI32ChPAT2_T1w.nii.gz", "_ses-bl_acq-ADNI32Ch_T1w.nii.gz",
                     "_ses-bl_acq-ADNI_T1w.nii.gz", "_ses-bl2_acq-ADNI32ChPAT2_T1w.nii.gz"]:
            if os.path.isfile(dir_bl + "anat/" + participant + suffix):
                  t1_bl = dir_bl + "anat/" + participant + suffix
      flair_bl = ""
      for suffix in ["_ses-bl2_FLAIR.nii.gz", "_ses-bl_FLAIR.nii.gz"]:
            if os.path.isfile(dir_bl + "anat/" + participant + suffix):
                  flair_bl = dir_bl + "anat/" + participant + suffix
      # there are also some naming inconsistencies at FU that have to be accounted for
      t1_fu = ""
      for suffix in ["_ses-fu_acq-ADNI32ChPAT2_T1w.nii.gz", "_ses-fu2_acq-ADNI32ChPAT2_T1w.nii.gz"]:
            if os.path.isfile(dir_fu + "anat/" + participant + suffix):
                  t1_fu = dir_fu + "anat/" + participant + suffix
      flair_fu = ""
      for suffix in ["_ses-fu2_FLAIR.nii.gz", "_ses-fu_FLAIR.nii.gz"]:
            if os.path.isfile(dir_fu + "anat/" + participant + suffix):
                  flair_fu = dir_fu + "anat/" + participant + suffix
      return dir_bl, dir_fu, t1_bl, t1_fu, flair_bl, flair_fu

# discern the participants with one timepoint that have both T1 and FLAIR and those with only FLAIR/T1
# they need different preprocessing steps
onetpboth = [] # n = 1583
onetpt1 = [] # n = 6
onetpflair = [] # n = 0
rest1tp = [] # n = 7 - these participants have neither T1 nor FLAIR scan
# we must also discern the participants with two timepoints that have 1 / 2 modalities at baseline / at follow-up
# 9 different combinations could exist
twotpbothboth = [] # all modalities at all timepoints --- n = 1068
twotpbotht1 = [] # all modalities at BL and only T1 at FU --- n = 1 --- only one subject - will thus not be processed via slurm
# SAMSEG expects that the same combination of contrasts is provided for each time point, so this subject will not be processed via the longitudinal pipeline
twotpbothflair = [] # all modalities at BL and only FLAIR at FU ---  n = 0
twotpt1both = [] # only t1 at BL and all modalities at FU --- n = 1 --- only one subject - will thus not be processed via slurm
# SAMSEG expects that the same combination of contrasts is provided for each time point, so this subject will not be processed via the longitudinal pipeline
twotpt1t1 = [] # only t1 at both timepoints ---  n = 0
twotpt1flair = [] # only t1 at baseline and only flair at FU ---  n = 0
twotpflairboth = [] # only flair at BL and all modalities at FU ---  n = 0
twotpflairt1 = [] # only flair at bl and only t1 at fu ---  n = 0
twotpflairflair = [] # only flair at all timepoints ---  n = 0
rest2tps = [] # n = 0
for participant, sessions in participants.items():
      dir_bl, dir_fu, t1_bl, t1_fu, flair_bl, flair_fu = pathtovariable(participant, sessions)
      if "NA" in sessions: # first deal with participants that have only 1 tp
            if sessions[0] != "NA":
                        if os.path.exists(flair_bl) and os.path.isfile(t1_bl):
                              onetpboth.append(participant)
                        elif os.path.isfile(flair_bl):
                              onetpflair.append(participant)
                        elif os.path.isfile(t1_bl):
                              onetpt1.append(participant)
                        else:
                              rest1tp.append(participant)
                  elif sessions[1] != "NA":
                        if os.path.isfile(flair_fu) and os.path.isfile(t1_fu):
                              onetpboth.append(participant)
                        elif os.path.isfile(flair_fu):
                              onetpflair.append(participant)
                        elif os.path.isfile(t1_fu):
                              onetpt1.append(participant)
                        else:
                              rest1tp.append(participant)
      else: # now handle participants with two timepoints
            if os.path.exists(t1_bl):
                  if os.path.exists(flair_bl):
                        if os.path.exists(t1_fu):
                              if os.path.exists(flair_fu):
                                    twotpbothboth.append(participant)
                              else:
                                    twotpbotht1.append(participant)
                        elif os.path.exists(flair_fu):
                              twotpbothflair.append(participant)
                        else:
                              rest2tps.append(participant)
                  elif os.path.exists(t1_fu):
                        if os.path.exists(flair_fu):
                              twotpt1both.append(participant)
                        else:
                              twotpt1t1.append(participant)
                  elif os.path.exists(flair_fu):
                        twotpt1flair.append(participant)
                  else:
                        rest2tps.append(participant)
            elif os.path.exists(flair_bl):
                  if os.path.exists(t1_fu):
                        if os.path.exists(flair_fu):
                              twotpflairboth.append(participant)
                        else:
                              twotpflairt1.append(participant)
                  elif os.path.exists(flair_fu):
                        twotpflairflair.append(participant)
                  else:
                        rest2tps.append(participant)
            else:
                  rest2tps.append(participant)

# save variables and functions for later use
dill.dump_session("/data/pt_life_freesurfer/samseg/scripts/prepared_vars.pkl")

# lastly we will create the directories that we will need
onetp = onetpt1 + onetpboth
for participant in onetp:
      if participants[participant][0] != "NA":
            if not os.path.exists("/data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/"):
                  os.mkdir("/data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/")
      else:
            if not os.path.exists("/data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/"):
                  os.mkdir("/data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/")

twotp = twotpbothboth + twotpt1both + twotpbotht1
for participant in twotp:
      if not os.path.exists("/data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/"):
                  os.mkdir("/data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/")
      if not os.path.exists("/data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/"):
                  os.mkdir("/data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/")
      if not os.path.exists("/data/pt_life_freesurfer/samseg/" + participant[4:] + "/"):
                  os.mkdir("/data/pt_life_freesurfer/samseg/" + participant[4:] + "/")
