# KaStat
Statisitical analysis for Kappa models. Currently, it carries out the following analyses.
1. Estimate the probability of a certain type of influence between two given Kappa rules.
2. Output the causal core that contributes most to the existence of the influence between the given Kappa rules.

# Compile 
----------------
mkdir build
cd build
cmake -DCMAKE_CXX_COMPILER=clang-omp++ ../src
make

# Usage
----------------
The command line is as follows:
        KaStat_sq <testfile> <modelfile> <ruleAname> <ruleBname> <influType> <KaSim> <simulationTime> <KaFlow>
        
where:
    - <testfile> is a text file containing a sequence of test specifications, give the path to it
    - <modelfile> is the file name and path of the Kappa model to be studied
    - <ruleAname> is the name of the rule that the user wants to study its impact on the other rule
    - <ruleBname> is the name of the rule that the user wants to study the impact from the other rule
    - <influType> is the type of the impact that the user wants to check, which can be either ``positive'' or ``negative''
    - <KaSim> is the KaSim executable, give the path to it
    - <simulationTime> is the limit of simulation time specified by the user
    - <KaFlow> is the KaFlow executable, give the path to it

