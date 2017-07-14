# OVERVIEW
This repository keeps the scripts and links to the data used to generate the results in SQUID manuscripts. There are 3 sections of results: on simulation data, on previously studied cell lines, on TCGA data. SQUID software is available at [https://github.com/Kingsford-Group/squid](https://github.com/Kingsford-Group/squid).

# ON SIMULATION DATA
### Software Used
Required:
- SQUID
- [STAR](https://github.com/alexdobin/STAR)
- [SpeedSeq](https://github.com/hall-lab/speedseq)
- [DELLY2](https://github.com/dellytools/delly)
- [LUMPY](https://github.com/arq5x/lumpy-sv)
- [Trinity](https://github.com/trinityrnaseq/trinityrnaseq)
- [Trans-ABySS](https://github.com/bcgsc/transabyss)
- [MUMMER3](https://github.com/marbl/MUMmer3)
- [Gmap](https://github.com/juliangehring/GMAP-GSNAP)

Optional (you will need them for simulating RNA-seq reads and SVs):
- [RSVsim](https://bioconductor.org/packages/release/bioc/html/RSVSim.html)
- [Flux Simulator](http://sammeth.net/confluence/display/SIM/Home)

The script build_simulation.sh will automatically download and install the above software.
```
./build_simulation.sh
```

Be sure to add these software into your path by running
```
export PATH=$(pwd)/tools/bin:${PATH}
```

### Workflow Description
![image of workflow with simulation data](doc/workflow_simulation.png)

### Re-producing the Result
- To skip the simulating step, download the simulated data simulation_data00 and simulation_data01 from [here](https://cmu.box.com/s/e9u6alp73rfdhfve2a51p6v391vweodq). Decompress them by running the following command in terminal.
```
cat simulation_data00 simulation_data01 > simulation_data.tar.gz
tar -xzvf simulation_data.tar.gz
```
To move on with the workflow, we provide a script to do alignment and TSV calling. You can generate all results by running
```
./script/runSimulationData.sh <decompressed data folder> <number of threads>
```

- If you want to simulate your own data, make sure the optional software are in your path. We write a script to run simulation as in the manuscript. Simply run the following command. (The whole simulation may take more than 1 day to finish, be patient ^ ^)
```
./script/runWholeSimulation.sh
```

Find the specification of output directory and files [here](doc/outputspec_simulation.md)

# ON PREVIOUSLY STUDIED CELL LINES
### Software Used

### Workflow Description

### Re-producing the Result

# ON TCGA DATA
### Software Used

### Workflow Description

### Re-producing the Result
