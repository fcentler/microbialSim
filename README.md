![µbialSim](microbialSimLogo300px.png)

# µbialSim

This is µbialSim (pronounced *microbialsim*), a dynamic Flux-Balance-Analysis-based simulator for complex microbial communities. Batch and chemostat operation can be simulated. Simulation results detail the temporal evolution of compound concentrations and biomass concentrations for each microbial species in the bioreactor, and internal flux distributions. For low diversity microbiomes Matlab's ODE solvers can be used, while for more complex communities an augmented forward Euler method is implemented in µbialSim.

## Citing

When using µbialSim, please cite [Popp, D. and F. Centler. 2020. “µbialSim: Constraint-Based Dynamic Simulation of Complex Microbiomes.”, Front. Bioeng. Biotechnol. 8 (June): 574. https://doi.org/10.3389/fbioe.2020.00574](https://www.frontiersin.org/articles/10.3389/fbioe.2020.00574/full).

## Getting Started

The latest version of µbialSim can be obtained via git:

```
git clone https://github.com/fcentler/microbialSim.git
```

### Prerequisites

µbialSim is implemented in Matlab and requires the COBRA Toolbox and/or CellNetAnalyzer for FBA computations. The examples make use of the former. As numerical solver, glpk, Gurobi and CPLEX can be used.

Information on the COBRA Toolbox is available here:

https://opencobra.github.io/cobratoolbox/stable/index.html

Information on CellNetAnalyzer is available here:

https://www2.mpi-magdeburg.mpg.de/projects/cna/cna.html

### Installing

For using the COBRA Toolbox, its location must be specified in the main simulator file `microbialSimMain.m` in lines 115ff, note that multiple locations can be specified, easing operation in different environments:

```
placesToLookForCobraToolbox = {
            'C:/cobratoolbox', ...
            'C:/Users/username/cobratoolbox', ...
            '/data/cobratoolbox' ...
            '../cobratoolbox' ...
            '/Users/username/cobratoolbox'
            };
```

If simulations are to be defined by external YAML files (see `config.yaml` for an example), [this Matlab yaml parser](https://github.com/MartinKoch123/yaml?tab=readme-ov-file) is required and must be copied into subdirectory `yaml` within the µbialSim directory. Simualtions can then be started as `microbialSimMain("config.yaml")`.

## Running the Examples

Three examples are included with µbialSim and can be run as follows.

### Example 1

The first example in which batch-culture growth of a single hydrogenotrophic species (*Methanococcus maripaludis*) is simulated can be run with the command:

```
microbialSimMain(1)
```

The FBA model is taken from *M. A. Richards, T. J. Lie, J. Zhang, S. W. Ragsdale, J. A. Leigh, and N. D. Price, “Exploring hydrogenotrophic methanogenesis: A genome scale metabolic reconstruction of Methanococcus maripaludis,” J. Bacteriol., vol. 198, no. 24, pp. 3379–3390, 2016.*

### Example 2

Batch-culture growth of a binary syntrophic, methanogenic community (*Syntrophobacter fumaroxidans* and *Methanospirillium hungatei*) transforming propionate to methane is started by:

```
microbialSimMain(2)
```

The FBA models are taken from *J. J. Hamilton, M. Calixto Contreras, and J. L. Reed, “Thermodynamics and H2 Transfer in a Methanogenic, Syntrophic Community,” PLOS Comput. Biol., vol. 11, no. 7, p. e1004364, 2015.*

### Example 3

Simulating batch growth of an eight species human gut microbiome (SIHUMIx consortium). Models are taken from the AGORA1 model collection containing 773 human gut microbiome species. Running this simulation first requires the unpacking of the file `AGORA-1.01-Western-Diet.zip` containing the AGORA models. The simulation can then be started by: 

```
microbialSimMain(3)
```

Any other subset of AGORA species, and even all 773 species together can be simulated. Models can either be selected by index in their alphabetical order, or by their filename. Note that for running the simulation with all 773 species, 64GB of RAM are necessary (loading models as a Matlab data structure after the initial loading as SBML files cuts memory demand in half). For complex microbiomes, consider deactivating the option to record limiting fluxes `recordLimitingFluxes` for faster simulation times. Speed can considerably be further improved by setting the parameter `maxDeviation` to `inf` in `microbialSimMain.m` at the expense of numerical accuracy.

The FBA models stem from *S. Magnúsdóttir et al., “Generation of genome-scale metabolic reconstructions for 773 members of the human gut microbiota,” Nat. Publ. Gr., vol. 35, no. 1, 2016.*

Models of the AGORA2 collection containing 7,302 models can be loaded as well. Individual models need to be [acquired from here](https://www.vmh.life/files/reconstructions/AGORA2/) and should be copied to the directory `models/AGORA-2.01-models`. Note that uptake constraints are unlimited in these models and require manual curation.

For chemostat simulations, the compound mix of the influx can be defined by a died file, as [obtainable from here](https://www.vmh.life/#nutrition). These files should be stored in the directory `models/diets`.


### Simulation Output

Three files are generated at the end of the simulation with a date and time stamp in the filename indicating the start of the simulation. All files hold Matlab data structures. The file `*_restartInit.mat` records the final state of the simulator and can be used as the initial conditions to continue the simulation in a subsequent run of µbialSim. The other file holds the simulated trajectory in the Matlab structure trajectory. The fields `time`, `compounds`, `biomass`, and `mu` hold the time, compound concentrations, biomass concentrations, and specific growth rates for each integration step. The field `FBA` stores data for each FBA model, including the temporal dynamics of all metabolic fluxes, and the mass balance for all exchange reactions. The field `limitingFluxes` stores over time those reactions, for which the current flux is identical to its current flux boundary (upper or lower). These are typically the growth-limiting fluxes. If only the evolution of biomass and compound concentrations are of interest, a simplified trajectory is saved to the file `*_dynamics.mat`. This file only records the evolution of concentrations and specific growth rates.
Additionally, once the simulation is finished, the trajectory is automatically visualized in Matlab figures detailing the evolution of compound and biomass concentrations as well as exchange fluxes for all species. To visualize only the most abundant microbial species or compounds, the function `plotTopBiomassCompounds()` can be used.

To inspect compound consumption and production rates at a given time, the function `getCmpndExchangeTable()` can be used, filtered (`filterCmpndExchangeTable()`), and saved to file as an adjacency matrix (`writeCmpndExchangeTableToFile()`) for subsequent visualization as a network.

## Known Issues

Parallel execution occasionally (but reproducible) hangs when running on AMD Ryzen CPUs.

## Built With

* [Matlab](https://www.mathworks.com/products/matlab.html) - Implementation language and numerical ODE solver
* [COBRA Toolbox](https://opencobra.github.io/cobratoolbox/stable/index.html) - For loading SBML models and running FBA computations
* [CellNetAnalyzer](https://www2.mpi-magdeburg.mpg.de/projects/cna/cna.html) - For loading SBML models and running FBA computations

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/fcentler/microbialSim/tags). Versions prior to v1.1.1 available at [https://git.ufz.de/UMBSysBio/microbialsim/](https://git.ufz.de/UMBSysBio/microbialsim/).

## Version History

* v1.2.0 - March 2024 - support for loading AGORA2 models, simulation parameters can now alternatively be set in an external YAML file, substrate concentration in inflow for chemostat simulations can be defined by external file, simplified trajectory files, support for numerical solvers Gurobi and CPLEX added, improved flexibility when restarting simulations, introduction of log levels 
* v1.1.1 - April 2020 - new plot routine for visualizing the most abundant microbial species and compounds only, new parameter for the numerical solver, minor bugfixes
* v1.1.0 - January 2020 - improved parallel performance, new function to estimate maximal uptake flux required to achieve given specific growth rate, new functions to collect, filter, and write compound exchange within the community at a given simulation time point, new feature to record growth-limiting compounds for each species during the simulation, bugfixes
* v1.0.0 - July 2019 - initial version

## Authors

* **Florian Centler**, University of Siegen, Faculty of Life Sciences, Siegen, Germany
* **Denny Popp**, UFZ - Helmholtz Centre for Environmental Research, Leipzig, Germany

Contact: florian.centler@uni-siegen.de

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details

## Acknowledgments

* This code was made possible by funding from the German Federal Ministry of Education and Research to Florian Centler (e:Bio project McBiogas, FKZ 031A317)
* Logo design by Rodrigo A. Colpo, UFZ - Helmholtz Centre for Environmental Research, Leipzig, Germany
