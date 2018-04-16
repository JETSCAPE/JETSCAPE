import re
import glob
import os

pattern_1 = r' [a-zA-Z0-9]+\_\w*\('
pattern_2 = r' [a-zA-Z0-9]+\_\w*\('
pattern_3 = r' [a-zA-Z0-9]+\_\w*\('

dirs = os.listdir('./')

excluded_files = ['AdSCFT.h', 'Matter.h', 'LBT.h', 'Martini.h', 'FourVector.h', 'PartonShower.h', 'JetClass.h', 'JetScapeParticles.h', 'HydroFromFile.h']

files = []
for fold in dirs:
    if not os.path.isdir(fold): continue
    all_files = os.listdir(fold)
    for f in all_files:
        if '.h' in f:
            # only read function names from header files
            files.append(os.path.join(fold, f))


func_names_before_change = []
for f in files:
    with open(f, 'r') as fin:
        # skip the folders in exclude_dirs
        fname = f.split('/')[1]
        if fname in excluded_files: continue
        fstr = fin.read()
        entries = re.findall(pattern_1, fstr) 
        entries += re.findall(pattern_2, fstr)
        entries += re.findall(pattern_3, fstr)
        if len(entries) != 0:
            print("\n=============================\n")
            print("For file: %s \n"%f)
            for entry in entries:
                print(entry[:-1]+'()',)
                func_name = entry[1:-1]
                '''skip some function names that from std/external libraries'''
                if func_name == 'make_unique': continue
                if func_name == 'sample_nucleons': continue
                if func_name == 'clean_hydro_event': continue
                func_names_before_change.append(func_name)

func_names_before_change = set(func_names_before_change)
print('func_names_before_change=', func_names_before_change)


files_to_change = []
for fold in (dirs + ['../examples/', '../jail/']):
    if not os.path.isdir(fold): continue
    all_files = os.listdir(fold)
    for f in all_files:
        if '.h' in f or '.cc' in f or '.cxx' in f:
            # only read function names from header files
            files_to_change.append(os.path.join(fold, f))
print(files_to_change)


def new_func_name_for(old_name):
    new_name = old_name
    if new_name[-1] == '_':
        new_name = new_name[:-1]
    #new_name = new_name.replace('_', ' ')
    return new_name.title().replace('_', '')

replacing_maps = {}
for name in func_names_before_change:
    replacing_maps[name] = new_func_name_for(name)

print(replacing_maps)

for f in files_to_change:
    f_new_str = ''
    with open(f, 'r') as fin:
        # skip the folders in exclude_dirs
        fstr = fin.read()
        f_new_str = fstr
        for old_str, new_str in replacing_maps.items():
            f_new_str = f_new_str.replace(old_str, new_str)

    with open(f, 'w') as fout:
        fout.write(f_new_str)

