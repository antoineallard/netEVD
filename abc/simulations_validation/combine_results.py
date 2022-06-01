import os
import glob
import pickle

dirs = glob.glob("../rough_results/simulations_validation_*")


file = "explored_parameters.txt"
for dir in dirs:
     if os.path.isfile(dir + "/" + file):
         # print(dir + "/" + file)
         os.remove(dir + "/" + file)


file = "simulation_results.txt"
cmd = "head -n 1 " + dirs[0] + "/" + file + " > " + "results/" + file
# print(cmd)
os.system(cmd)
cmd = ["awk 'FNR>1'"]
for dir in dirs:
    cmd.append(dir + "/" + file)
cmd.append(">>")
cmd.append("results/" + file)
# print(" ".join(cmd))
os.system(" ".join(cmd))


file = "refined_simulation_results.txt"
# print(cmd)
cmd = ["awk 'FNR>1'"]
for dir in dirs:
    if os.path.isfile(dir + "/" + file):
        cmd2 = "head -n 1 " + dir + "/" + file + " > " + "results/" + file
        os.system(cmd2)
        cmd.append(dir + "/" + file)
cmd.append(">>")
cmd.append("results/" + file)
# print(" ".join(cmd))
os.system(" ".join(cmd))


with open(dirs[0] + '/config_dict.pickle', 'rb') as file:
    configuration = pickle.load(file)
for dir in dirs[1:]:
    with open(dir + '/config_dict.pickle', 'rb') as file:
        meta_data = pickle.load(file)
        configuration["N_simulations"] += meta_data["N_simulations"]
configuration['path'] =  'results'
configuration['data_path'] =  "validation_data/time_series.txt"
configuration['folder_path'] = ''
configuration['folder_name'] = ''
with open("results/" + 'config_dict.pickle', 'wb') as file:
    pickle.dump(configuration, file, protocol=pickle.HIGHEST_PROTOCOL)
