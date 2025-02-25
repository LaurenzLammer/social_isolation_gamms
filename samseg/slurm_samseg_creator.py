# this script will create multiple slurm scripts to run freesurfer's samseg across all life subjects
# to be run parallelized via the slurm compute cluster
# if participants have both a T1 scan and a FLAIR scan both will be used
# if participants underwent MRIs at baseline and follow-up, the longitudinal version of samseg will be run
# in the longitudinal version a template based on all scans has to be created before the individual timepoints can be analysed
import os
import dill

# load the data prepared with slurm_samseg_creator_preparation.py
dill.load_session("/data/pt_life_freesurfer/samseg/scripts/prepared_vars.pkl")



# for participants with both modalities, their FALIRS have to be coregistered and reformatted to the T1 scans' format
# this is achieved by running two consecutive scripts
# prepare a slurm script for n participants with one timepoint and two modalities for co-registration
n = 583 # how many participants should be processed via the script
m = 1000 # starting with which participant from the 1tp list
c = 10 # how many cores should be requested
# create header
slurmheader(c = c, n = n, m = m, name = "slurm_onetp_coregistration_")

for participant in onetpboth[m:n+m]:
      dir_bl, dir_fu, t1_bl, t1_fu, flair_bl, flair_fu = pathtovariable(participant, participants[participant])
      if participants[participant][0] != "NA":
            print("srun --exclusive -n 1 --output /data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/coregistration_log.txt " +
                  "freesurfer --version 7.4.1 mri_coreg --mov " + flair_bl + " --ref " + t1_bl + " --reg /data/pt_life_freesurfer/samseg/" +
                  participant[4:] + "_bl/flairToT1.lta &",
                  file=open(("/data/pt_life_freesurfer/samseg/scripts/slurm_onetp_coregistration_" + str(m) + "_to_" + str(m+n-1)), "a"))
      else:
            print("srun --exclusive -n 1 --output /data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/coregistration_log.txt " +
                  "freesurfer --version 7.4.1 mri_coreg --mov " + flair_fu + " --ref " + t1_fu + " --reg /data/pt_life_freesurfer/samseg/" +
                  participant[4:] + "_fu/flairToT1.lta &",
                  file=open(("/data/pt_life_freesurfer/samseg/scripts/slurm_onetp_coregistration_" + str(m) + "_to_" + str(m+n-1)), "a"))
# add wait at end of the commands
print("wait", file=open(("/data/pt_life_freesurfer/samseg/scripts/slurm_onetp_coregistration_" + str(m) + "_to_" + str(m+n-1)), "a"))

# we will now verify that the coregistrations have not failed
for participant in onetpboth[1000:]:
      if participants[participant][0] != "NA":
            with open('/data/pt_life_freesurfer/samseg/' + participant[4:] + "_bl/coregistration_log.txt") as f:
                  lines = [line for line in f]
      else:
            with open('/data/pt_life_freesurfer/samseg/' + participant[4:] + "_fu/coregistration_log.txt") as f:
                  lines = [line for line in f]
      if "mri_coreg done" not in lines[-2]:
            print(participant)
            print(participant, file=open(("/data/pt_life_freesurfer/samseg/scripts/onetp_coregistration_errors" + str(
                        m) + "_to_" + str(m + n - 1)), "a"))


# prepare a slurm script for n participants with one timepoint and two modalities for reformatting
n = 583 # how many participants should be processed via the script
m = 1000 # starting with which participant from the 1tp list
c = 40 # how many cores should be requested
# create header
slurmheader(c = c, n = n, m = m, name = "slurm_onetp_reformatting_")

for participant in onetpboth[m:n+m]:
      dir_bl, dir_fu, t1_bl, t1_fu, flair_bl, flair_fu = pathtovariable(participant, participants[participant])
      if participants[participant][0] != "NA":
            print("srun --exclusive -n 1 --output /data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/reformatting_log.txt " +
                  "freesurfer --version 7.4.1 mri_vol2vol --mov " + flair_bl + " --reg /data/pt_life_freesurfer/samseg/" +
                  participant[4:] + "_bl/flairToT1.lta " + " --targ " + t1_bl + " --o /data/pt_life_freesurfer/samseg/" +
                  participant[4:] + "_bl/flair_reg.nii &",
                  file=open(("/data/pt_life_freesurfer/samseg/scripts/slurm_onetp_reformatting_" + str(m) + "_to_" + str(m+n-1)), "a"))
      else:
            print("srun --exclusive -n 1 --output /data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/reformatting_log.txt " +
                  "freesurfer --version 7.4.1 mri_vol2vol --mov " + flair_fu + " --reg /data/pt_life_freesurfer/samseg/" +
                  participant[4:] + "_fu/flairToT1.lta " + " --targ " + t1_fu + " --o /data/pt_life_freesurfer/samseg/" +
                  participant[4:] + "_fu/flair_reg.nii &",
                  file=open(("/data/pt_life_freesurfer/samseg/scripts/slurm_onetp_reformatting_" + str(m) + "_to_" + str(m+n-1)), "a"))
# add wait at end of the commands
print("wait", file=open(("/data/pt_life_freesurfer/samseg/scripts/slurm_onetp_reformatting_" + str(m) + "_to_" + str(m+n-1)), "a"))

# we will now verify that the reformattings have not failed
for participant in onetpboth[1000:]:
      if participants[participant][0] != "NA":
            with open('/data/pt_life_freesurfer/samseg/' + participant[4:] + "_bl/reformatting_log.txt") as f:
                  lines = [line for line in f]
      else:
            with open('/data/pt_life_freesurfer/samseg/' + participant[4:] + "_fu/reformatting_log.txt") as f:
                  lines = [line for line in f]
      if "mri_vol2vol done" not in lines[-1]:
            print(participant)
            print(participant, file=open(("/data/pt_life_freesurfer/samseg/scripts/onetp_reformatting_errors" + str(
                        m) + "_to_" + str(m + n - 1)), "a"))

# once their FLAIRS are coregistered and reformatted, the participants with one tp and both modalities can be processed by samseg
# prepare a slurm script for n participants with one timepoint and both modalities to run the samseg algorithm
n = 283 # how many participants should be processed via the script
m = 1300 # starting with which participant from the 1tp list
c = 40 # how many cores should be requested
# create header
slurmheader(c = c, n = n, m = m, name = "slurm_cross_1tp_both")
## add the commands for the individual participants
for participant in onetpboth[m:m+n]:
      dir_bl, dir_fu, t1_bl, t1_fu, flair_bl, flair_fu = pathtovariable(participant, participants[participant])
      # first participants for which only baseline scans exists
      if participants[participant][0] != "NA":
            print("srun --exclusive -n 1 --output /data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/log.txt " +
                  "freesurfer --version 7.4.1 run_samseg --input " + t1_bl +
                  " /data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/flair_reg.nii" +
                  " --pallidum-separate --lesion --lesion-mask-pattern 0 1 --output /data/pt_life_freesurfer/samseg/" +
                        participant[4:] + "_bl/ --threads 1 &",
                        file=open(("/data/pt_life_freesurfer/samseg/scripts/" + "slurm_cross_1tp_both" + str(m) + "_to_" + str(m+n-1)), "a"))
      else:
            print("srun --exclusive -n 1 --output /data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/log.txt " +
                  "freesurfer --version 7.4.1 run_samseg --input " + t1_fu +
                  " " + " /data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/flair_reg.nii" +
                  " --pallidum-separate --lesion --lesion-mask-pattern 0 1 --output /data/pt_life_freesurfer/samseg/" +
                  participant[4:] + "_fu/ --threads 1 &",
                  file=open(("/data/pt_life_freesurfer/samseg/scripts/" + "slurm_cross_1tp_both" + str(m) + "_to_" + str(m + n - 1)),
                            "a"))
# add wait at end of the commands
print("wait", file=open(("/data/pt_life_freesurfer/samseg/scripts/" + "slurm_cross_1tp_both" + str(m) + "_to_" + str(m + n - 1)), "a"))

# we will now verify that the samseg algorithms have not failed
for participant in onetpboth:
      if participants[participant][0] != "NA":
            with open('/data/pt_life_freesurfer/samseg/' + participant[4:] + "_bl/log.txt") as f:
                  lines = [line for line in f]
      else:
            with open('/data/pt_life_freesurfer/samseg/' + participant[4:] + "_fu/log.txt") as f:
                  lines = [line for line in f]
      if "run_samseg complete" not in lines[-1]:
            print(participant)
            print(participant, file=open(("/data/pt_life_freesurfer/samseg/scripts/onetp_samseg_both_errors" + str(
                        m) + "_to_" + str(m + n - 1)), "a"))

# prepare a slurm script for n participants with one timepoint and only T1 to run the samseg algorithm
n = 6 # how many participants should be processed via the script
m = 0 # starting with which participant from the 1tp list
c = 6 # how many cores should be requested
# create header
slurmheader(c = c, n = n, m = m, name = "slurm_cross_1tp_t1")
## add the commands for the individual participants
for participant in onetpt1[m:m+n]:
      dir_bl, dir_fu, t1_bl, t1_fu, flair_bl, flair_fu = pathtovariable(participant, participants[participant])
      # first participants for which only baseline scans exists
      if participants[participant][0] != "NA":
            print("srun --exclusive -n 1 --output /data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/log.txt " +
                        "freesurfer --version 7.4.1 run_samseg --input " + t1_bl +
                        " --output /data/pt_life_freesurfer/samseg/" +
                        participant[4:] + "_bl/ --threads 1 &",
                        file=open(("/data/pt_life_freesurfer/samseg/scripts/" + "slurm_cross_1tp_t1" + str(m) + "_to_" + str(m+n-1)), "a"))
      else:
            print("srun --exclusive -n 1 --output /data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/log.txt " +
                  "freesurfer --version 7.4.1 run_samseg --input " + t1_fu +
                  " --output /data/pt_life_freesurfer/samseg/" +
                  participant[4:] + "_fu/ --threads 1 &",
                  file=open(("/data/pt_life_freesurfer/samseg/scripts/" + "slurm_cross_1tp_t1" + str(m) + "_to_" + str(m + n - 1)),
                            "a"))
# add wait at end of the commands
print("wait", file=open(("/data/pt_life_freesurfer/samseg/scripts/" + "slurm_cross_1tp_t1" + str(m) + "_to_" + str(m + n - 1)), "a"))

# we will now verify that the samseg algorithms have not failed
for participant in onetpt1:
      if participants[participant][0] != "NA":
            with open('/data/pt_life_freesurfer/samseg/' + participant[4:] + "_bl/log.txt") as f:
                  lines = [line for line in f]
      else:
            with open('/data/pt_life_freesurfer/samseg/' + participant[4:] + "_fu/log.txt") as f:
                  lines = [line for line in f]
      if "run_samseg complete" not in lines[-1]:
            print(participant)
            print(participant, file=open(("/data/pt_life_freesurfer/samseg/scripts/onetp_samseg_t1_errors" + str(
                        m) + "_to_" + str(m + n - 1)), "a"))

# for participants with measurements at two timepoints, templates must be created to make use of the longitudinal pipeline
# this must take place before flairs are coregistered and reformatted

# create a script for template creation for participants with two timepoints and both modalities at each timepoint
n = 68 # how many participants should be processed via the script
m = 1000 # starting with which participant from the 1tp list
c = 10 # how many cores should be requested
# create header
slurmheader(c = c, n = n, m = m, name = "slurm_twotp_template_creation_")
for participant in twotpbothboth[m:n+m]:
      dir_bl, dir_fu, t1_bl, t1_fu, flair_bl, flair_fu = pathtovariable(participant, participants[participant])
      print("srun --exclusive -n 1 --output /data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/template_creation_log.txt " +
            "freesurfer --version 7.4.1 mri_robust_template --mov " + t1_bl + " " + t1_fu +  " --template /data/pt_life_freesurfer/samseg/" +
            participant[4:] + "_bl/mean.mgz --satit --mapmov " + "/data/pt_life_freesurfer/samseg/" + participant[4:] +
            "_bl/t1_reg.mgz /data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/t1_reg.mgz",
            file=open(("/data/pt_life_freesurfer/samseg/scripts/slurm_twotp_template_creation_" + str(m) + "_to_" +
                       str(m + n - 1)), "a"))

# we will now verify that the template creations have not failed
for participant in twotpbothboth[200:]:
      with open("/data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/template_creation_log.txt") as f:
                  lines = [line for line in f]
      if "registration took" not in lines[-10]:
            print(participant)
            print(participant, file=open(("/data/pt_life_freesurfer/samseg/scripts/twotp_template_errors" + str(
                        m) + "_to_" + str(m + n - 1)), "a"))

# create a script for coregistration for participants with two timepoints and both modalities at each timepoint
n = 868 # how many participants should be processed via the script
m = 200 # starting with which participant from the 1tp list
c = 10 # how many cores should be requested - here c may be 2 * n because the flairs of both timepoints have to be coregistered
# create header
slurmheader(c = c, n = n, m = m, name = "slurm_twotp_coregistration_")
for participant in twotpbothboth[m:n+m]:
      dir_bl, dir_fu, t1_bl, t1_fu, flair_bl, flair_fu = pathtovariable(participant, participants[participant])
      print("srun --exclusive -n 1 --output /data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/coregistration_log.txt " +
            "freesurfer --version 7.4.1 mri_coreg --mov " + flair_bl + " --ref /data/pt_life_freesurfer/samseg/" + participant[4:] +
            "_bl/t1_reg.mgz --reg /data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/flairToT1.lta &",
            file=open(("/data/pt_life_freesurfer/samseg/scripts/slurm_twotp_coregistration_" + str(m) + "_to_" +
            str(m + n - 1)), "a"))
      print("srun --exclusive -n 1 --output /data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/coregistration_log.txt " +
            "freesurfer --version 7.4.1 mri_coreg --mov " + flair_fu + " --ref /data/pt_life_freesurfer/samseg/" + participant[4:] +
            "_fu/t1_reg.mgz --reg /data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/flairToT1.lta &",
            file=open(("/data/pt_life_freesurfer/samseg/scripts/slurm_twotp_coregistration_" + str(m) + "_to_" +
                       str(m + n - 1)), "a"))
# add wait at end of the commands
print("wait", file=open(("/data/pt_life_freesurfer/samseg/scripts/slurm_twotp_coregistration_" + str(m) + "_to_" + str(m + n - 1)), "a"))

# we will now verify that the coregistrations have not failed
for participant in twotpbothboth[200:]:
      with open("/data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/coregistration_log.txt") as f:
                  lines = [line for line in f]
      if "mri_coreg done" not in lines[-2]:
            print(participant + "_bl", file=open(("/data/pt_life_freesurfer/samseg/scripts/twotp_coregistration_errors" + str(
                        m) + "_to_" + str(m + n - 1)), "a"))
      with open("/data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/coregistration_log.txt") as f:
                  lines = [line for line in f]
      if "mri_coreg done" not in lines[-2]:
            print(participant)
            print(participant + "_fu", file=open(("/data/pt_life_freesurfer/samseg/scripts/twotp_coregistration_errors" + str(
                        m) + "_to_" + str(m + n - 1)), "a"))

# create a script for reformatting for participants with two timepoints and both modalities at each timepoint
n = 868 # how many participants should be processed via the script
m = 200 # starting with which participant from the 1tp list
c = 50 # how many cores should be requested - here c = 2 * n because the flairs of both timepoints have to be coregistered
# create header
slurmheader(c = c, n = n, m = m, name = "slurm_twotp_reformatting_")
for participant in twotpbothboth[m:n+m]:
      dir_bl, dir_fu, t1_bl, t1_fu, flair_bl, flair_fu = pathtovariable(participant, participants[participant])
      print("srun --exclusive -n 1 --output /data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/reformatting_log.txt " +
            "freesurfer --version 7.4.1 mri_vol2vol --mov " + flair_bl + " --reg /data/pt_life_freesurfer/samseg/" +
            participant[4:] + "_bl/flairToT1.lta " + " --targ /data/pt_life_freesurfer/samseg/" + participant[4:] +
            "_bl/t1_reg.mgz --o /data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/flair_reg.nii &",
            file=open(("/data/pt_life_freesurfer/samseg/scripts/slurm_twotp_reformatting_" + str(m) + "_to_" + str(m+n-1)), "a"))
      print("srun --exclusive -n 1 --output /data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/reformatting_log.txt " +
            "freesurfer --version 7.4.1 mri_vol2vol --mov " + flair_fu + " --reg /data/pt_life_freesurfer/samseg/" +
            participant[4:] + "_fu/flairToT1.lta " + " --targ /data/pt_life_freesurfer/samseg/" + participant[4:] +
            "_fu/t1_reg.mgz --o /data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/flair_reg.nii &",
            file=open(("/data/pt_life_freesurfer/samseg/scripts/slurm_twotp_reformatting_" + str(m) + "_to_" + str(
                  m + n - 1)), "a"))
# add wait at end of the commands
print("wait", file=open(("/data/pt_life_freesurfer/samseg/scripts/slurm_twotp_reformatting_" + str(m) + "_to_" + str(m + n - 1)), "a"))

# we will now verify that the reformattings have not failed
for participant in twotpbothboth[200:]:
      with open("/data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/reformatting_log.txt") as f:
                  lines = [line for line in f]
      if "mri_vol2vol done" not in lines[-1]:
            print(participant)
            #print(participant + "_bl", file=open(("/data/pt_life_freesurfer/samseg/scripts/twotp_reformatting_errors" + str(
             #           m) + "_to_" + str(m + n - 1)), "a"))
      with open("/data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/reformatting_log.txt") as f:
                  lines = [line for line in f]
      if "mri_vol2vol done" not in lines[-1]:
            print(participant)
            print(participant + "_fu", file=open(("/data/pt_life_freesurfer/samseg/scripts/twotp_reformatting_errors" + str(
                        m) + "_to_" + str(m + n - 1)), "a"))

# after the template has been created and the flair scans have been coregistered and reformatted
# the longitudinal samseg algorithm can be run

n = 283 # how many participants should be processed via the script
m = 500 # starting with which participant from the 1tp list
c = 100 # how many cores should be requested
# create header
slurmheader(c = c, n = n, m = m, name = "slurm_twotp_samseg_long_", mem = "20")
for participant in twotpbothboth[m:n+m]:
      dir_bl, dir_fu, t1_bl, t1_fu, flair_bl, flair_fu = pathtovariable(participant, participants[participant])
      print("srun --exclusive -n 1 --output /data/pt_life_freesurfer/samseg/" + participant[4:] + "/samseg_long_log.txt " +
            "freesurfer --version 7.4.1 run_samseg_long --timepoint /data/pt_life_freesurfer/samseg/" + participant[4:] +
            "_bl/t1_reg.mgz /data/pt_life_freesurfer/samseg/" + participant[4:] + "_bl/flair_reg.nii --timepoint /data/pt_life_freesurfer/samseg/" +
            participant[4:] + "_fu/t1_reg.mgz /data/pt_life_freesurfer/samseg/" + participant[4:] + "_fu/flair_reg.nii" +
            " --pallidum-separate --lesion --lesion-mask-pattern 0 1 --output /data/pt_life_freesurfer/samseg/" +
            participant[4:] + "/ --threads 1 &",
            file=open(("/data/pt_life_freesurfer/samseg/scripts/slurm_twotp_samseg_long_" + str(m) + "_to_" + str(m + n - 1)), "a"))
# add wait at end of the commands
print("wait", file=open(("/data/pt_life_freesurfer/samseg/scripts/slurm_twotp_samseg_long_" + str(m) + "_to_" + str(m + n - 1)), "a"))

# we will now verify that the longitudinal samseg algorithms have not failed
todo = []
for participant in twotpbothboth:
      if os.path.exists("/data/pt_life_freesurfer/samseg/" + participant[4:] + "/samseg_long_log.txt"):
            n += 1
      with open("/data/pt_life_freesurfer/samseg/" + participant[4:] + "/samseg_long_log.txt") as f:
                  lines = [line for line in f]
      if "run_samseg_long complete" not in lines[-1]:
            print(participant)
            todo.append(participant)
            print(participant, file=open(("/data/pt_life_freesurfer/samseg/scripts/twotp_samseg_long_errors" + str(m) +
                                          "_to_" + str(m + n - 1)), "a")
# there is only one participant in twotpbotht1 and twotpt1both each so I will not run them via slurm
# due to the varying number of available modalities the longitudinal pipeline can't be run
# each timepoint will be treated like an individual participant when it comes to Samseg


