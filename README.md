OptStoic python package
========================
Perform optStoic analysis using Python code that share the same data files with GAMS code.

Note: All the examples are specific for glycolysis pathway generation. 

## Install
- Next, setup a virtual environment in Python 2.7.
```bash
# Install the virtualenv package
pip install virtualenv
# Test virtualenv
virtualenv --version
# Create a project folder
cd project_folder
# Create a virtual environment call optstoic_env
virtualenv optstoic_env
# Activate your environment
source optstoic_env/bin/activate
```

- Then, install one of the solvers in the following section.


- Next, clone this repository in your `project_folder` and setup. This should install all the Python dependencies.
```
cd project_folder
git clone https://github.com/maranasgroup/optstoic-python.git
cd optstoic-python
python setup.py install
```

- To run nosetests after setup, you need to run setup as below:
```
python setup.py test
```

## Requirement
At least one of the following optimization solvers should be installed. To solve the loopless optStoic formulation, an optimization solver other than GLPK is recommended.

1. GLPK 4.47 installation
    ```bash
    wget  http://ftp.gnu.org/gnu/glpk/glpk-4.47.tar.gz
    tar -xvzf glpk-4.47.tar.gz
    cd  ~/glpk-4.47
    ./configure
    make
    make install
    #if the program is successfully installed, you should get an output by typing
    glpsol --version
    ```

2. GUROBI Optimization provide academic license for free (https://www.gurobi.com/). Install gurobipy following the instruction provided by GUROBI. 

3. [SCIP Optimization Suite](https://scip.zib.de/) >= v4.0.0
```
sudo apt-get install libgmp-dev libreadline-dev zlib1g-dev libncurses5-dev
tar xvf scipoptsuite-6.0.0.tgz
cd scipoptsuite-6.0.0/
make
make test
cd scip-6.0.0/
sudo make install INSTALLDIR="/usr/local/"
/usr/local/bin/scip --version
```

4. [CPLEX Optimizer](https://www.ibm.com/analytics/cplex-optimizer)

## Current project dependencies
1. [PuLP](https://github.com/coin-or/pulp). Run the [test](https://www.coin-or.org/PuLP/main/installing_pulp_at_home.html#testing-your-pulp-installation).

2. Graphviz (Optional, for drawing pathway). The [Graphviz](https://www.graphviz.org/) software is required before installing the graphviz python package. 
    ```bash
    #If you have root access
    sudo apt-get install graphviz

    #If you do not have root access (you can get a different version of Graphviz from their website https://www.graphviz.org/download/)
    cd $HOME
    mkdir -p bin/graphviz
    wget http://www.graphviz.org/pub/graphviz/stable/SOURCES/graphviz-2.38.0.tar.gz
    tar xvf graphviz-2.38.0.tar.gz
    cd graphviz-2.38.0
    ./configure --prefix=$HOME/bin/graphviz
    make && make install
    # Check if the graphviz is working
    cd $HOME/bin/graphviz/bin
    dot -V
    # Add the following line to your .bashrc
    export PATH=$PATH:$HOME/bin/graphviz/bin

    #Install the Python graphviz package
    pip install graphviz
    ```

3. [Component-Contribution](https://github.com/eladnoor/component-contribution) (*Optional, unless you want to perform MDF analysis)

## Tests
After cloning the repo or setup, please run tests as followed:
```
pip install nose 
nosetests -s -v 
```

To test the optimization algorithms, please run tests as followed. The runtime depends on the solvers selected by PuLP.
```
from optstoicpy.test.testAll import test_all_optimization_scripts
test_all_optimization_scripts()
```

## Usage
Read the [tutorial](https://github.com/maranasgroup/optstoic-python/blob/master/optstoicpy/examples/methods.md).

## Development
To continue development with the code, please create a virtual environment and use `python setup.py develop` for installation.

## Reference
Please cite [Ng, C.Y., Wang, L., Chowdhury, A. et al. Pareto Optimality Explanation of the Glycolytic Alternatives in Nature. Sci Rep 9, 2633 (2019). https://doi.org/10.1038/s41598-019-38836-9](https://www.nature.com/articles/s41598-019-38836-9).
