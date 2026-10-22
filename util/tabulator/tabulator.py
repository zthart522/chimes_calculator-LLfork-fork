import sys
import os
import glob

sys.path.append(os.path.normpath(os.getcwd()))

if not os.path.exists("config.py"):
    
    print("Error: Cannot find pes_generator config file")
    print("       Searched for: config.py in the current directory")
    exit(0)

# Read the config file

import config
    
os.system( "python3" + " " + config.CHMS_REPO + "util/pes_generator/src/pes_generator.py")

generated_files = glob.glob("chimes_scan_2b.type_*dat")

for i in range(len(generated_files)):
    
    infile        = generated_files[i]
    output_energy = generated_files[i] + ".energy"
    output_force  = generated_files[i] + ".force"
    
    # Count the number of lines 
    
    with open(infile, 'r') as ifstream:
        line_count = 0
        for line in ifstream:
            line_count += 1
    
    # Write the out the tabulated files
    
    energystream = open(output_energy,'w')
    forcestream  = open(output_force,'w')
    
    energystream.write(str(line_count-1) + '\n')
    forcestream.write(str(line_count-1) + '\n')
    
    first = True

    with open(infile, 'r') as ifstream:
        for line in ifstream:
            if (first):
                first = False
                continue


           # energystream.write(line.split()[0] + " 0" + '\n')
            energystream.write(line.split()[0] + " " + line.split()[1] + '\n')
           # forcestream .write(line.split()[0] + " 0" + '\n')
            forcestream .write(line.split()[0] + " " + line.split()[2] + '\n')

    
    energystream.close()
    forcestream .close()

### Trial 3b generation ### 

generated_files = glob.glob("chimes_scan_3b.type_*dat")

for i in range(len(generated_files)):
    
    infile        = generated_files[i]
    output_energy = generated_files[i] + ".energy"
    output_force  = generated_files[i] + ".force"
    
    # Count the number of lines 
    
    with open(infile, 'r') as ifstream:
        line_count = 0
        for line in ifstream:
            line_count += 1
    
    # Write the out the tabulated files
    
    energystream = open(output_energy,'w')
    forcestream  = open(output_force,'w')
    
    energystream.write(str(line_count-1) + '\n')
    forcestream.write(str(line_count-1) + '\n')
    
    first = True

    with open(infile, 'r') as ifstream:
        for line in ifstream:
            if (first):
                first = False
                continue


            energystream.write(line.split()[0] + " " + line.split()[1]+ " " + line.split()[2] + " " + line.split()[3] + '\n')
           # energystream.write(line.split()[0] + " " + line.split()[1]+ " " + line.split()[2] + " 0" + '\n')
            forcestream .write(line.split()[0] + " " + line.split()[1]+ " " + line.split()[2] + " " + line.split()[4]+ " " + line.split()[5] + " " + line.split()[6] + '\n')
           # forcestream .write(line.split()[0] + " " + line.split()[1]+ " " + line.split()[2] + " 0 0 0"  + '\n')

    
    energystream.close()
    forcestream .close()

# Write new parameter file if name is specified
if len(sys.argv) == 2:
    
    # Check that params are not being overwritten, then write to files
    if os.path.realpath(config.PARAM_FILE) == os.path.realpath(sys.argv[1]):
        print("\nERROR: Overriding parameter file in config.py is not allowed.")
        print("Rerun script with different tabulated parameter filename or paste lines manually.\n")
        sys.exit(1)

    # Open original params and tabulated version
    paramstream = open(config.PARAM_FILE,'r')
    tabstream   = open(sys.argv[1],'w')
        
    # Adjust 2B lines (PAIRTYPE PARAMS: {index} {E1} {E2})
    generated_files = glob.glob("chimes_scan_2b.type_*dat")
    for i in range(len(generated_files)):
        infile        = generated_files[i]
        index         = infile.split("type_")[1].split(".dat")[0]

        for line in paramstream:
            if line.strip().split(":")[0] == "PAIRTYPE PARAMS" and index == line.strip().split()[2]:
                tabstream.write(line.rstrip("\n") + f" TABULATED /{infile}\n")
                break
            else:
                tabstream.write(line)
            
    # Adjust 3B lines (INDEX: {index} ATOMS: {E1} {E2} {E3})
    generated_files = glob.glob("chimes_scan_3b.type_*dat")
    for i in range(len(generated_files)):
        infile        = generated_files[i]
        index         = infile.split("type_")[1].split(".dat")[0]

        for line in paramstream:
            print(line.strip().split(":")[0])
            if line.strip().split(":")[0] == "INDEX" and "ATOMS" in line.strip().split(":")[1] and index == line.strip().split()[1]:
                tabstream.write(line.rstrip("\n") + f" TABULATED /{infile}\n")
                break
            else:
                tabstream.write(line)

    # Finish writing to the end of the file
    for line in paramstream:
        tabstream.write(line)

    paramstream.close()
    tabstream.close()
