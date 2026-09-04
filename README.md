This repository contains the scripts needed to perform the main analyses
described in the associated paper ("qbaconfound: A flexible Monte Carlo
probabilistic bias analysis to unmeasured confounding" submitted to Statistical
Methods in Medical Research). It includes scripts to simulate datasets and then
analyse these datasets for the simulation study ("Simulation study" folder
corresponding to Section 3 in the paper), and scripts to download and clean
data, and analyse it for the applied example ("Applied example" folder
corresponding to Section 4 in the paper). It also contains "Bayesian function"
which contains a modified function of the *unm_glm* function of R package
*unmconf*.

The structure of the repository follows as:

- Simulation study
  - Data simulation
    - Generates datasets for scenarios A and B of simulation study I.do:
      generates 500 simulated datasets for scenarios A and B of simulation study I.
      Datasets are saved in a folder called "Data" on pathway "Simulation study
      I/Scenario A/" or "Simulation study I/Scenario B/"
    - Generates datasets for scenario G of simulation study IV.do: generates 500
      simulated datasets for scenario G of simulation study IV. Datasets are
      saved in a folder called "Data" on pathway "Simulation study IV/yst_xcon/"
    - Scenario A - simulates a dataset.do: simulates a single dataset for
      scenario A of simulation studies I and II
    - Scenario B - simulates a dataset.do: simulates a single dataset for
      scenario B of simulation studies I and II
    - Scenario C - simulates a dataset.do: simulates a single dataset for
      scenario C of simulation study III
    - Scenario D - simulates a dataset.do: simulates a single dataset for
      scenario D of simulation study III
    - Scenario E - simulates a dataset.do: simulates a single dataset for
      scenario E of simulation study IV
    - Scenario F - simulates a dataset.do: simulates a single dataset for
      scenario F of simulation study IV

  - Fits MCQBA
    - Scenario A - fits MCQBA.do: applies a Monte Carlo QBA to data with a
      binary outcome and continuous unmeasured confounder (scenario A)
    - Scenario B - fits MCQBA.do: applies a Monte Carlo QBA to data with a
      continuous and two continuous unmeasured confounders (scenario B)
    - Scenario C - fits MCQBA.do: applies a Monte Carlo QBA to data with a
      binary outcome and binary unmeasured confounder (scenario C)
    - Scenario D - fits MCQBA.do: applies a Monte Carlo QBA to data with a
      continuous outcome and two binary unmeasured confounder (scenario D)
    - Scenario E - fits MCQBA.do: applies a Monte Carlo QBA to data with a
      continuous outcome and binary and continuous unmeasured confounders
      (scenario E)
    - Scenario F - fits MCQBA.do: applies a Monte Carlo QBA to data with a
      nominal outcome and continuous unmeasured confounder (scenario F)
    - Scenario G - fits MCQBA.do: applies a Monte Carlo QBA to data with a
      survival outcome and continuous unmeasured confounder (scenario G)

  - Simulation study I
    - Scenario A
      - Data
        - 500 CSV files from Dataset_1.csv to Dataset_500.csv of the simulated
          datasets for scenario A of simulation study I
      - Monte Carlo
        - Do files
          - Scenario A - logit - runs sim study I for Monte Carlo QBA.do: runs
            simulation study I for Monte Carlo QBA for scenario A
        - Results
          - Empty folder to store the results from running the above do file
      - Bayesian
        - Seeds
          - Generates random integers.do: generates 8,000 random integers for
            seeds
          - RandomIntegers.csv: lists 8,000 random integers
        - R files
          - ScenarioA-glm-RunsSimStudyIforBayesianQBA.R: runs simulation study I
            for Bayesian QBA for scenario A
        - Results
          - Empty folder to store the results from running the above R file
    - Scenario B
      - Data
        - 500 CSV files from Dataset_1.csv to Dataset_500.csv of the simulated
          datasets for scenario B of simulation study I
      - Monte Carlo
        - Do files
          - Scenario B - regress - runs sim study I for Monte Carlo QBA.do: runs
            simulation study I for Monte Carlo QBA for scenario B
        - Results
          - Empty folder to store the results from running the above do file
      - Bayesian
        - Seeds
          - Generates random integers.do: generates 8,000 random integers for
            seeds
          - RandomIntegers.csv: lists 8,000 random integers
        - R files
          - ScenarioB-lm-RunsSimStudyIforBayesianQBA.R: runs simulation study I
            for Bayesian QBA for scenario B
        - Results
          - Empty folder to store the results from running the above R file

  - Simulation study II
    - Do files
      - Scenario A - logit - runs sim study II for Monte Carlo QBA.do: runs
        simulation study II for Bayesian QBA for scenario A
      - Scenario B - regress - runs sim study II for Monte Carlo QBA.do: runs
        simulation study II for Bayesian QBA for scenario B
    - Results
      - Empty folder to store the results from running the above do files

  - Simulation study III
    - Do files
      - Scenario C - logit - runs sim study III for Monte Carlo QBA.do: runs
        simulation study III for Bayesian QBA for scenario C
      - Scenario D - regress - runs sim study III for Monte Carlo QBA.do: runs
        simulation study III for Bayesian QBA for scenario D
    - Results
      - Empty folder to store the results from running the above do files

  - Simulation study IV
    - yst_xcon
      - 500 CSV files from Dataset_1.csv to Dataset_500.csv of the simulated
        datasets for scenario G of simulation study IV
    - Do files
      - Scenario E - regress - runs sim study IV for Monte Carlo QBA.do: runs
        simulation study IV for Bayesian QBA for scenario E
      - Scenario F - mlogit - runs sim study IV for Monte Carlo QBA.do: runs
        simulation study IV for Bayesian QBA for scenario F
      - Scenario G - Cox PH - runs sim study IV for Monte Carlo QBA.do: runs
        simulation study IV for Bayesian QBA for scenario G
    - Results
      - Empty folder to store the results from running the above do files

- Applied example
  - GeneratesDataForNHANESExample.R: generates data for the NHANES example
  - Data
    - Data.csv: cleaned NHANES data for analysis
  - Monte Carlo
    - Apply MCQBA to NHANES example using prior for piM.do: runs Monte Carlo QBA
      with a prior distribution for the marginal prevalence
    - Apply MCQBA to NHANES example using prior for the intercept.do: runs Monte
      Carlo QBA with a prior distribution for the intercept of logistic
      regression
  - Bayesian
    - ApplyBayesianQBAToNHANESExample.R: runs Bayesian QBA

- Bayesian function
  - MyUnm_glm.R: modified *unm_glm* function from R package *unmconf* 1.0.0 for
    running Bayesian QBA
