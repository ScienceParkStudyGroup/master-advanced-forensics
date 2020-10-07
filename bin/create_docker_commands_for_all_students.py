#!/usr/bin/env python3


# Usage: create_docker_commands_for_all_students.py [file with student to virtual machine correspondence]

# Input file should have tabulated separated values (.tsv or .txt)
# Input file format should contain these columns with this naming scheme.
# student	     machine	  password	port
# Maura Cook 	 machine-01	  maura	    8787
# Reini van Hal  machine-02	  reini	    8788
# ...

import pandas as pd
import sys

students_to_machine = sys.argv[1]

df = pd.read_csv(student_to_machine, sep="\t")

def create_docker_command(row):
    """Takes a Pandas row and return the docker command with corresponding student name + pwd + port number"""	
    student =     row["student"]
    machine_nb =  row["machine"]
    password   =  row["password"]
    port =        row["port"]

    docker_cmd = "docker run --name " + machine_nb + " -e PASSWORD=" + password + " -p 8787:" + str(port) + " scienceparkstudygroup/master-gls:openr-latest"
    print(docker_cmd)
    return docker_cmd

docker_commands = df.apply(create_docker_command, axis = 1, result_type='reduce')

