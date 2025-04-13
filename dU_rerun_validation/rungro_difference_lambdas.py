import subprocess
import numpy as np

gromacs_bin = ${location of double precision gromacs binary}
# Ex: gromacs_bin = "/Users/lwhitm1/projects/paper/pkgs/gromacs_double/bin/gmx_d"

# these can be any up to 10, as that is the number of states defined in the file.
initial_lambda_state = 0
final_lambda_states = [0,1,2,3,4,5,6,7,8,9]  # always include 0!!
nsteps = 5 # set from the MDP
system_top = "reference.top"
system_gro = "reference.gro"
template_mdp = "md_template.mdp"
#dispcorr = "EnerPres" # options  = "EnerPres" and "No"
dispcorr = "No" # options  = "EnerPres" and "No"
#vdw_modifier = "Potential-switch" # options = "Potential-switch", "Potential-shift", and "Force-switch"
#vdw_modifier = "Force-switch" # options = "Potential-switch", "Potential-shift", and "Force-switch"
#vdw_modifier = "Potential-shift" # options = "Potential-switch", "Potential-shift", and "Force-switch"
vdw_modifier = "Potential-switch" # options = "Potential-switch", "Potential-shift", and "Force-switch"

if dispcorr == "No":
    potential_index = 10
elif dispcorr == "EnerPres":
    potential_index = 11


def replace_lines(template,fe_state,output_file):

    f = open(template,'r')
    mdp_template_lines = f.readlines()
    f.close()

    f = open(output_file,'w')
    for line in mdp_template_lines:
        if line[0] == ';':
            continue;
        if "REPLACE_INIT_LAMBDA_STATE" in line:
            print(f"init-lambda-state = {fe_state:d}",file=f)
        elif "REPLACE_DISPCORR" in line:
            print(f"dispcorr = {dispcorr}",file=f)
        elif "REPLACE_VDW_MODIFIER" in line:
            print(f"vdw-modifier = {vdw_modifier}",file=f)
        elif "vdw-switch" in line:
            if vdw_modifier != "potential-switch":
                print(line,file=f)
        else:
            print(line,file=f)
    f.close()

def read_dhdl(dhdl_file, states):

    du_vals = list()
    f = open(dhdl_file,'r')
    lines = f.readlines()
    f.close()
    for line in lines:
        if line[0] == '#' or line[0] == "@":
            continue
        vals = line.split()
        du_vals_sub = list()
        for i in states:
            du_vals_sub.append(float(vals[i+4]))
        du_vals.append(du_vals_sub)
    return np.array(du_vals)

def read_edr(edr_file, index):

    edr_vals = list()
    f = open(edr_file,'r')
    lines = f.readlines()
    f.close()
    for line in lines:
        if line[0] == '#' or line[0] == "@":
            continue
        vals = line.split()
        edr_vals.append(float(vals[index+1]))
    return np.array(edr_vals)
    
# first, do the initial run.  Create the new mdp.
prefix = "fe_run" + str(initial_lambda_state)
replace_lines("md_template.mdp", initial_lambda_state, prefix+".mdp")
subprocess.run([gromacs_bin,"grompp","-f",prefix+".mdp", "-p", system_top, "-c", system_gro, "-o", prefix+".tpr"])
subprocess.run([gromacs_bin,"mdrun","-nt","1","-s",prefix+".tpr", "-g", prefix+".log", "-e",  prefix + ".edr", "-o", prefix+"out.trr", "-c", prefix + "out.gro", "-dhdl", prefix + ".dhdl.xvg"])
# read the .edr, since it's higher precision
# pipe the quantities of interest here
selection_process = subprocess.Popen(["echo", str(potential_index), "0"], stdout=subprocess.PIPE)
subprocess.run([gromacs_bin,"energy", "-dp", "-f", prefix + ".edr", "-o", prefix + ".energy.xvg"], stdin = selection_process.stdout)

val_init_pot = read_edr(prefix + ".energy.xvg", index = 0)

# now, read the edr file and extract the potential energy differences of interest
dvals_run = read_dhdl(prefix + ".dhdl.xvg", states = final_lambda_states)

# this is an array of length "number of frames output x number of final lambdas x number of frames"

# save old prefix
prefix_run = prefix
vals_eval = np.zeros([nsteps,len(final_lambda_states)])
for ip, i in enumerate(final_lambda_states):
    #next, run the rerun, using the trajectory from above
    prefix = "fe_eval" + str(i)
    replace_lines(template_mdp, i, prefix+".mdp")
    subprocess.run([gromacs_bin,"grompp","-f",prefix+".mdp", "-p", system_top, "-c", system_gro, "-o", prefix+".tpr"])
    subprocess.run([gromacs_bin,"mdrun","-nt","1","-rerun", prefix_run + "out.trr", "-s",prefix+".tpr", "-g", prefix+".log", "-e",  prefix + ".edr", "-o", prefix+"out.trr", "-c", prefix + "out.gro", "-dhdl", prefix + ".dhdl.xvg"])
    # pipe the quantities of interest here

    # read the .edr, since it's higher precision
    # pipe the quantities of interest here
    selection_process = subprocess.Popen(["echo", str(potential_index), "0"], stdout=subprocess.PIPE)
    subprocess.run([gromacs_bin, "energy", "-dp", "-f", prefix + ".edr", "-o", prefix + ".energy.xvg"], stdin = selection_process.stdout)    
    vals_eval[:,ip] = read_edr(prefix + ".energy.xvg", index = 0)

dvals_eval = (vals_eval.transpose() - val_init_pot).transpose()

#What we find is that with rerun, there is always a configuration
#dependent offset in each frame.  But it is NOT lambda dependent, will
#not affect the free energies. So we need to correct this by
#subtracting off the lambda = 0 potential for each configuration.

corr_dvals_eval = (dvals_eval.transpose() - dvals_eval[:,0]).transpose()


# now output results.
print("run")
for j in range(len(final_lambda_states)):
    print(f"{j:14d}",end="")
print("")
for i in range(nsteps):
    for j in range(len(final_lambda_states)):
        print(f"{dvals_run[i,j]:14.8f}",end="")
    print("")

print("evaluated")
for j in range(len(final_lambda_states)):
    print(f"{j:14d}",end="")
print("")
for i in range(nsteps):
    for j in range(len(final_lambda_states)):
        print(f"{corr_dvals_eval[i,j]:14.8f}",end="")
    print("")

print("difference")
for j in range(len(final_lambda_states)):
    print(f"{j:14d}",end="")
print("")
for i in range(nsteps):
    for j in range(len(final_lambda_states)):
        print(f"{dvals_run[i,j]-corr_dvals_eval[i,j]:14.8f}",end="")
    print("")

    
