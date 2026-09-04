This repository contains the scripts needed to perform the main analyses
described in the following paper:

**A flexible Monte Carlo quantitative bias analysis for unmeasured confounding**

It includes scripts to simulate datasets and then analyse these datasets for the
simulation study ([SimulationStudy/](./SimulationStudy/) corresponding to
Section 3 in the paper), and scripts to download and clean data, and analyse it
for the applied example ([AppliedExample/](./AppliedExample/) corresponding to
Section 4 in the paper). It also contains [Bayesian
function/](./BayesianFunction/) which contains a modified function of the
*unm_glm* function of R package *unmconf*.

The structure of the repository follows as:

- [**Simulationtudy/**](./SimulationStudy/)

  - [Master do file for running Monte Carlo QBA in simulation studies I to
    IV.do](./SimulationStudy/Master do file for running Monte Carlo QBA in simulation studies I to IV.do):
    Master do file for running simulation studies I to IV for Monte Carlo QBA.
    Runs do files from folder "Do files" and post results to a folder called
    "Results".

  - [**Data simulation/**](./SimulationStudy/Data simulation/)
    - [Generates datasets for scenarios A and B of simulation study
      I.do](./SimulationStudy/Data simulation/Generates datasets for scenarios A and B of simulation study I.do):
      Generates 500 simulated datasets for scenarios A and B of simulation study
      I and saves in [Simulation study I/Scenario
      A/Data/](./SimulationStudy/Simulation study I/Scenario A/Data/) and
      [Simulation study I/Scenario
      B/Data/](./SimulationStudy/Simulation study I/Scenario B/Data/)
      respectively.

    - [Generates datasets for scenario G of simulation study
      IV.do](./SimulationStudy/Data simulation/Generates datasets for scenario G of simulation study IV.do):
      Generates 500 simulated datasets for scenario G of simulation study IV and
      saves in [Simulation study
      IV/Data/](./SimulationStudy/Simulation study IV/Data/).

    - [Scenario A - simulates a
      dataset.do](./SimulationStudy/Data simulation/Scenario
      A - simulates a dataset.do):
      Simulates a single dataset for scenario A of simulation studies I, II and VI.

    - [Scenario B - simulates a
      dataset.do](./SimulationStudy/Data simulation/Scenario
      B - simulates a dataset.do):
      Simulates a single dataset for scenario B of simulation studies I, II and VI.

    - [Scenario C - simulates a
      dataset.do](./SimulationStudy/Data simulation/Scenario
      C - simulates a dataset.do):
      Simulates a single dataset for scenario C of simulation study III.

    - [Scenario D - simulates a
      dataset.do](./SimulationStudy/Data simulation/Scenario
      D - simulates a dataset.do):
      Simulates a single dataset for scenario D of simulation study III.

    - [Scenario E - simulates a
      dataset.do](./SimulationStudy/Data simulation/Scenario
      E - simulates a dataset.do):
      Simulates a single dataset for scenario E of simulation study IV.

    - [Scenario F - simulates a
      dataset.do](./SimulationStudy/Data simulation/Scenario F - simulates a dataset.do):
      Simulates a single dataset for scenario F of simulation study IV.

  - [**Fits MCQBA/**](./SimulationStudy/Fits MCQBA/)
    - [Scenario A - fits
      MCQBA.do](./SimulationStudy/Fits MCQBA/Scenario A - fits MCQBA.do):
      Applies Monte Carlo QBA to data with a binary outcome and continuous
      unmeasured confounder (scenario A).

    - [Scenario B - fits
      MCQBA.do](./SimulationStudy/Fits MCQBA/Scenario
      B - fits MCQBA.do):
      Applies Monte Carlo QBA to data with a continuous and two continuous
      unmeasured confounders (scenario B).

    - [Scenario C - fits
      MCQBA.do](./SimulationStudy/Fits MCQBA/Scenario
      C - fits MCQBA.do):
      Applies Monte Carlo QBA to data with a binary outcome and binary
      unmeasured confounder (scenario C).

    - [Scenario D - fits
      MCQBA.do](./SimulationStudy/Fits MCQBA/Scenario
      D - fits MCQBA.do):
      Applies Monte Carlo QBA to data with a binary outcome and two binary
      unmeasured confounder (scenario D).

    - [Scenario E - fits
      MCQBA.do](./SimulationStudy/Fits MCQBA/Scenario
      E - fits MCQBA.do):
      Applies Monte Carlo QBA to data with a continuous outcome and two
      unmeasured confounders (one binary and one continuous) (scenario E).

    - [Scenario F - fits
      MCQBA.do](./SimulationStudy/Fits MCQBA/Scenario
      F - fits MCQBA.do):
      Applies Monte Carlo QBA to data with a nominal outcome and continuous
      unmeasured confounder (scenario F).

    - [Scenario G - fits
      MCQBA.do](./SimulationStudy/Fits MCQBA/Scenario G - fits MCQBA.do):
      Applies Monte Carlo QBA to data with a survival outcome and continuous
      unmeasured confounder (scenario G).

  - [**Simulation study I/**](./SimulationStudy/Simulation study I/)
    - [**Scenario A/**](./SimulationStudy/Simulation study I/Scenario A)
      - [**Data/**](./SimulationStudy/Simulation study I/Scenario A/Data/)
        - 500 csv files from
          [Dataset_1.csv](./SimulationStudy/Simulation study I/Scenario A/Data/Dataset_1.csv)
          to
          [Dataset_500.csv](./SimulationStudy/Simulation study I/Scenario A/Data/Dataset_500.csv)
          of the simulated datasets for scenario A of simulation studies I, II
          and VI
      - [**Monte
        Carlo/**](./SimulationStudy/Simulation study I/Scenario A/Monte Carlo/)
        - [**Do
          files/**](./SimulationStudy/Simulation study I/Scenario A/Monte Carlo/Do files/)
          - [Scenario A - logit - runs sim study I for Monte Carlo
            QBA.do](./SimulationStudy/Simulation study I/Scenario A/Monte Carlo/Do files/Scenario A - logit - runs sim study I for Monte Carlo QBA.do):
            Runs scenario A of simulation study I for Monte Carlo QBA.
        - [**Results/**](./SimulationStudy/Simulation study I/Scenario A/Monte Carlo/Results/)
          - empty folder to store the results from running the above do file
      - [**Bayesian/**](./SimulationStudy/Simulation study I/Scenario A/Bayesian/)
        - [**Seeds/**](./SimulationStudy/Simulation study I/Scenario A/Bayesian/Seeds/)
          - [Generates random
            integers.do](./SimulationStudy/Simulation study I/Scenario A/Bayesian/Seeds/Generates random integers.do):
            Generates 8,000 random integers for seeds.
          - [RandomIntegers.csv](./SimulationStudy/Simulation study I/Scenario A/Bayesian/Seeds/RandomIntegers.csv):
            Lists 8,000 random integers.
        - [**R
          files/**](./SimulationStudy/Simulation study I/Scenario A/Bayesian/R files/)
          - [ScenarioA-glm-RunsSimStudyIforBayesianQBA.R](./SimulationStudy/Simulation study I/Scenario A/Bayesian/R files/ScenarioA-glm-RunsSimStudyIforBayesianQBA.R) :
            Runs scenario A of simulation study I for Bayesian QBA.
        - [**Results/**](./SimulationStudy/Simulation study I/Scenario A/Bayesian/Results/)
          - empty folder to store the results from running the above R file
    - [**Scenario B/**](./SimulationStudy/Simulation study I/Scenario B/)
      - [**Data/**](./SimulationStudy/Simulation study I/Scenario B/Data/)
        - 500 csv files from
          [Dataset_1.csv](./SimulationStudy/Simulation study I/Scenario B/Data/Dataset_1.csv)
          to
          [Dataset_500.csv](./SimulationStudy/Simulation study I/Scenario B/Data/Dataset_500.csv)
          of the simulated datasets for scenario B of simulation studies I, II
          and VI
      - [**Monte
        Carlo/**](./SimulationStudy/Simulation study I/Scenario B/Monte Carlo)
        - [**Do
          files/**](./SimulationStudy/Simulation study I/Scenario B/Monte Carlo/Do files/)
          - [Scenario B - regress - runs sim study I for Monte Carlo
            QBA.do](./SimulationStudy/Simulation study I/Scenario B/Monte Carlo/Do files/Scenario B - regress - runs sim study I for Monte Carlo QBA.do):
            Runs scenario B of simulation study I for Monte Carlo QBA.
        - [**Results/**](./SimulationStudy/Simulation study I/Scenario B/Monte Carlo/Results/)
          - empty folder to store the results from running the above do file
      - [**Bayesian/**](./SimulationStudy/Simulation study I/Scenario B/Bayesian)
        - [**Seeds/**](./SimulationStudy/Simulation study I/Scenario B/Bayesian/Seeds/)
          - [Generates random
            integers.do](./SimulationStudy/Simulation study I/Scenario B/Bayesian/Seeds/Generates random integers.do):
            Generates 8,000 random integers for seeds.
          - [RandomIntegers.csv](./SimulationStudy/Simulation study I/Scenario B/Bayesian/Seeds/RandomIntegers.csv):
            Lists 8,000 random integers.
        - [**R
          files/**](./SimulationStudy/Simulation study I/Scenario B/Bayesian/R files/)
          - [ScenarioB-lm-RunsSimStudyIforBayesianQBA.R](./SimulationStudy/Simulation study I/Scenario B/Bayesian/R files/ScenarioB-lm-RunsSimStudyIforBayesianQBA.R):
            Runs scenario B of simulation study I for Bayesian QBA.
        - [**Results/**](./SimulationStudy/Simulation study I/Scenario B/Bayesian/Results/)
          - empty folder to store the results from running the above R file

  - [**Simulation study II/**](./SimulationStudy/Simulation study II/)
    - [**Do files/**](./SimulationStudy/Simulation study II/Do files/)
      - [Scenario A - logit - runs sim study II for Monte Carlo
        QBA.do](./SimulationStudy/Simulation study II/Do files/Scenario A - logit - runs sim study II for Monte Carlo QBA.do):
        Runs scenario A of simulation study II for Monte Carlo QBA.
      - [Scenario B - regress - runs sim study II for Monte Carlo
        QBA.do](./SimulationStudy/Simulation study II/Do files/Scenario B - regress - runs sim study II for Monte Carlo QBA.do):
        Runs scenario B of simulation study II for Monte Carlo QBA.
    - [**Results/**](./SimulationStudy/Simulation study II/Results/)
      - empty folder to store the results from running the above do files

  - [**Simulation study III/**](./SimulationStudy/Simulation study III/)
    - [**Do files/**](./SimulationStudy/Simulation study III/Do files/)
      - [Scenario C - logit - runs sim study III for Monte Carlo
        QBA.do](./SimulationStudy/Simulation study III/Do files/Scenario C - logit - runs sim study III for Monte Carlo QBA.do):
        Runs scenario C of simulation study III for Monte Carlo QBA.
      - [Scenario D - regress - runs sim study III for Monte Carlo
        QBA.do](./SimulationStudy/Simulation study III/Do files/Scenario D - regress - runs sim study III for Monte Carlo QBA.do):
        Runs scenario D of simulation study III for Monte Carlo QBA.
    - [**Results/**](./SimulationStudy/Simulation study III/Results/)
      - empty folder to store the results from running the above do files

  - [**Simulation study IV/**](./SimulationStudy/Simulation study IV/)
    - [**Data/**](./SimulationStudy/Simulation study IV/Data/)
      - 500 csv files from
        [Dataset_1.csv](./SimulationStudy/Simulation study IV/Data/Dataset_1.csv)
        to
        [Dataset_500.csv](./SimulationStudy/Simulation study IV/Data/Dataset_500.csv)
        of the simulated datasets for scenario G of simulation study IV
    - [**Do files/**](./SimulationStudy/Simulation study IV/Do files/)
      - [Scenario E - regress - runs sim study IV for Monte Carlo
        QBA.do](./SimulationStudy/Simulation study IV/Do files/Scenario E - regress - runs sim study IV for Monte Carlo QBA.do):
        Runs scenario E of simulation study IV for Monte Carlo QBA.
      - [Scenario F - mlogit - runs sim study IV for Monte Carlo
        QBA.do](./SimulationStudy/Simulation study IV/Do files/Scenario F - mlogit - runs sim study IV for Monte Carlo QBA.do):
        Runs scenario F of simulation study IV for Monte Carlo QBA.
      - [Scenario G - Cox PH - runs sim study IV for Monte Carlo
        QBA.do](./SimulationStudy/Simulation study IV/Do files/Scenario G - Cox PH - runs sim study IV for Monte Carlo QBA.do):
        Runs scenario G of simulation study IV for Monte Carlo QBA.
    - [**Results/**](./SimulationStudy/Simulation study IV/Results/)
      - empty folder to store the results from running the above do files

  - [**Simulation study V/**](./SimulationStudy/Simulation study V/)
    - [**Data**](./SimulationStudy/Simulation study V/Data/)
      - 500 dta files from
        [Dataset_1.dta](./SimulationStudy/Simulation study V/Data/Dataset_1.dta)
        to
        [Dataset_500.dta](./SimulationStudy/Simulation study V/Data/Dataset_500.dta)
        of the simulated datasets for scenarios A and B of simulation study V
    - [**Do files/**](./SimulationStudy/Simulation study V/Do files/)
      - [Scenario A - logit - runs sim study V for Monte Carlo
        QBA.do](./SimulationStudy/Simulation study V/Do files/Scenario A - logit - runs sim study V for Monte Carlo QBA.do):
        Runs scenario A of simulation study V for Monte Carlo QBA.
      - [Scenario B - regress - runs sim study V for Monte Carlo
        QBA.do](./SimulationStudy/Simulation study V/Do files/Scenario B - regress - runs sim study V for Monte Carlo QBA.do):
        Runs scenario B of simulation study V for Monte Carlo QBA.
    - [**Results/**](./SimulationStudy/Simulation study V/Results/)
      - empty folder to store the results from running the above do files

  - [**Simulation study VI/**](./SimulationStudy/Simulation study VI/)
    - [**Do files/**](./SimulationStudy/Simulation study VI/Do files/)
      - [Scenario A - logit - runs sim study VI for Monte Carlo
        QBA.do](./SimulationStudy/Simulation study VI/Do files/Scenario A - logit - runs sim study VI for Monte Carlo QBA.do):
        Runs scenario A of simulation study VI for Monte Carlo QBA.
      - [Scenario B - regress - runs sim study VI for Monte Carlo
        QBA.do](./SimulationStudy/Simulation study VI/Do files/Scenario B - regress - runs sim study VI for Monte Carlo QBA.do):
        Runs scenario B of simulation study VI for Monte Carlo QBA.
    - [**Results/**](./SimulationStudy/Simulation study VI/Results/)
      - empty folder to store the results from running the above do files

- [**AppliedExample/**](./AppliedExample/)
  - [GeneratesDataForNHANESExample.R](./AppliedExample/GeneratesDataForNHANESExample.R):
    Generates data for the NHANES example.
  - [**Data/**](./AppliedExample/Data/)
    - [Data.csv](./AppliedExample/Data/Data.csv): cleaned NHANES data for
      analysis
  - [**Monte Carlo/**](./AppliedExample/Monte Carlo/)
    - [Apply MCQBA to NHANES example using prior for
      piM.do](./AppliedExample/Monte Carlo/Apply MCQBA to NHANES
      example using prior for piM.do):
      Applies Monte Carlo QBA with a prior distribution for the marginal
      prevalence.
    - [Apply MCQBA to NHANES example using prior for the
      intercept.do](./AppliedExample/Monte Carlo/Apply MCQBA to NHANES example using prior for the intercept.do):
      Applies Monte Carlo QBA with a prior distribution for the intercept of
      logistic regression.
    - [Apply MCQBA to NHANES example using inaccurate
      priors.do](./AppliedExample/Monte Carlo/Apply MCQBA to NHANES example using inaccurate priors.do):
      Applies Monte Carlo QBA using inaccurate priors.
  - [**Bayesian/**](./AppliedExample/Bayesian/)
    - [ApplyBayesianQBAToNHANESExample.R](./AppliedExample/Bayesian/ApplyBayesianQBAToNHANESExample.R):
      Applies Bayesian QBA.

- [**BayesianFunction/**](./BayesianFunction/)
  - [MyUnm_glm.R](./BayesianFunction/MyUnm_glm.R): modified *unm_glm* function
    from R package *unmconf* 1.0.0 for running Bayesian QBA

## Publications

Emily Kawabata, Chin Yang Shapland, Tom Palmer, David Carslake, Kate Tilling,
Rachael Hughes. A flexible Monte Carlo quantitative bias analysis for unmeasured
confounding. Statistical Methods in Medical Research (accepted in 2026).

See also:
[https://www.medrxiv.org/content/10.1101/2025.08.12.25333217v1](https://www.medrxiv.org/content/10.1101/2025.08.12.25333217v1)

## Further information

If you would like any further information, please contact
emily.kawabata@bristol.ac.uk.
