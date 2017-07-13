# KaStat

Statistical analysis for Kappa models. Currently, it carries out the following analyses.

1. Estimate the probability of a certain type of influence between two given Kappa rules.

2. Output the causal core that contributes most to the existence of the influence between the given Kappa rules.

# Compile 

    mkdir build
    cd build
    cmake -DCMAKE_CXX_COMPILER=clang-omp++ ../src
    make


# Usage

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

An example:
```
./KaStat_sq ../statistical_test/test01 ../models/egfr.ka "1/E-R" "2/R-R" positive ../KaSim 100 ../KaFlow
```

# Output

[given statistical test] [estimated prob] [# of trajectories holding the given type of influence] [# of total trajectories]
e.g.
BEST 0.01 0.99 1 1: estimate = 0.995633, successes = 227, samples = 227

If successes != 0, output:
"The story with this file name - " + maxratiofilename + " - contributes most to this influence with a ratio as " + maxratio + "%."

Otherwise, output:
"The statistical analysis indicates that there is no " + influType + " influence from " + string(argv[3]) [RuleA] + " to " + string(argv[4]) [RuleB].



