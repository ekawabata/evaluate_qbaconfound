This repository contains the scripts needed to perform the main analyses
described in the following paper:

**A flexible Monte Carlo quantitative bias analysis for unmeasured confounding**

It includes scripts to simulate datasets and then analyse these datasets for the
simulation study ([Simulation study/](./Simulation%20study/) corresponding to
Section 3 in the paper), and scripts to download and clean data, and analyse it
for the applied example ([Applied example/](./Applied%20example/) corresponding
to Section 4 in the paper). It also contains [Bayesian
function/](./Bayesian%20function/) which contains a modified function of the
*unm_glm* function of R package *unmconf*.

The structure of the repository follows as:

- [**Simulation study/**](./Simulation%20study/)

  - [Master do file for running Monte Carlo QBA in simulation studies I to
    IV.do](./Simulation%20study/Master%20do%20file%20for%20running%20Monte%20Carlo%20QBA%20in%20simulation%20studies%20I%20to%20IV.do):
    Master do file for running simulation studies I to IV for Monte Carlo QBA.
    Runs do files from folder "Do files" and post results to a folder called
    "Results".

  - [**Data simulation/**](./Simulation%20study/Data%20simulation/)
    - [Generates datasets for scenarios A and B of simulation study
      I.do](./Simulation%20study/Data%20simulation/Generates%20datasets%20for%20scenarios%20A%20and%20B%20of%20simulation%20study%20I.do):
      Generates 500 simulated datasets for scenarios A and B of simulation study
      I and saves in [Simulation study I/Scenario
      A/Data/](./Simulation%20study/Simulation%20study%20I/Scenario%20A/Data/)
      and [Simulation study I/Scenario
      B/Data/](./Simulation%20study/Simulation%20study%20I/Scenario%20B/Data/)
      respectively.

    - [Generates datasets for scenario G of simulation study
      IV.do](./Simulation%20study/Data%20simulation/Generates%20datasets%20for%20scenario%20G%20of%20simulation%20study%20IV.do):
      Generates 500 simulated datasets for scenario G of simulation study IV and
      saves in [Simulation study
      IV/Data/](./Simulation%20study/Simulation%20study%20IV/Data/).

    - [Scenario A - simulates a
      dataset.do](./Simulation%20study/Data%20simulation/Scenario%20A%20-%20simulates%20a%20dataset.do):
      Simulates a single dataset for scenario A of simulation studies I, II and VI.

    - [Scenario B - simulates a
      dataset.do](./Simulation%20study/Data%20simulation/Scenario%20B%20-%20simulates%20a%20dataset.do):
      Simulates a single dataset for scenario B of simulation studies I, II and VI.

    - [Scenario C - simulates a
      dataset.do](./Simulation%20study/Data%20simulation/Scenario%20C%20-%20simulates%20a%20dataset.do):
      Simulates a single dataset for scenario C of simulation study III.

    - [Scenario D - simulates a
      dataset.do](./Simulation%20study/Data%20simulation/Scenario%20D%20-%20simulates%20a%20dataset.do):
      Simulates a single dataset for scenario D of simulation study III.

    - [Scenario E - simulates a
      dataset.do](./Simulation%20study/Data%20simulation/Scenario%20E%20-%20simulates%20a%20dataset.do):
      Simulates a single dataset for scenario E of simulation study IV.

    - [Scenario F - simulates a
      dataset.do](./Simulation%20study/Data%20simulation/Scenario%20F%20-%20simulates%20a%20dataset.do):
      Simulates a single dataset for scenario F of simulation study IV.

  - [**Fits MCQBA/**](./Simulation%20study/Fits%20MCQBA/)
    - [Scenario A - fits
      MCQBA.do](./Simulation%20study/Fits%20MCQBA/Scenario%20A%20-%20fits%20MCQBA.do):
      Applies Monte Carlo QBA to data with a binary outcome and continuous
      unmeasured confounder (scenario A).

    - [Scenario B - fits
      MCQBA.do](./Simulation%20study/Fits%20MCQBA/Scenario%20B%20-%20fits%20MCQBA.do):
      Applies Monte Carlo QBA to data with a continuous and two continuous
      unmeasured confounders (scenario B).

    - [Scenario C - fits
      MCQBA.do](./Simulation%20study/Fits%20MCQBA/Scenario%20C%20-%20fits%20MCQBA.do):
      Applies Monte Carlo QBA to data with a binary outcome and binary
      unmeasured confounder (scenario C).

    - [Scenario D - fits
      MCQBA.do](./Simulation%20study/Fits%20MCQBA/Scenario%20D%20-%20fits%20MCQBA.do):
      Applies Monte Carlo QBA to data with a binary outcome and two binary
      unmeasured confounder (scenario D).

    - [Scenario E - fits
      MCQBA.do](./Simulation%20study/Fits%20MCQBA/Scenario%20E%20-%20fits%20MCQBA.do):
      Applies Monte Carlo QBA to data with a continuous outcome and two
      unmeasured confounders (one binary and one continuous) (scenario E).

    - [Scenario F - fits
      MCQBA.do](./Simulation%20study/Fits%20MCQBA/Scenario%20F%20-%20fits%20MCQBA.do):
      Applies Monte Carlo QBA to data with a nominal outcome and continuous
      unmeasured confounder (scenario F).

    - [Scenario G - fits
      MCQBA.do](./Simulation%20study/Fits%20MCQBA/Scenario%20G%20-%20fits%20MCQBA.do):
      Applies Monte Carlo QBA to data with a survival outcome and continuous
      unmeasured confounder (scenario G).

  - [**Simulation study I/**](./Simulation%20study/Simulation%20study%20I/)
    - [**Scenario
      A/**](./Simulation%20study/Simulation%20study%20I/Scenario%20A)
      - [**Data/**](./Simulation%20study/Simulation%20study%20I/Scenario%20A/Data/)
        - 500 csv files from
          [Dataset_1.csv](./Simulation%20study/Simulation%20study%20I/Scenario%20A/Data/Dataset_1.csv)
          to
          [Dataset_500.csv](./Simulation%20study/Simulation%20study%20I/Scenario%20A/Data/Dataset_500.csv)
          of the simulated datasets for scenario A of simulation studies I, II
          and VI
      - [**Monte
        Carlo/**](./Simulation%20study/Simulation%20study%20I/Scenario%20A/Monte%20Carlo/)
        - [**Do
          files/**](./Simulation%20study/Simulation%20study%20I/Scenario%20A/Monte%20Carlo/Do%20files/)
          - [Scenario A - logit - runs sim study I for Monte Carlo
            QBA.do](./Simulation%20study/Simulation%20study%20I/Scenario%20A/Monte%20Carlo/Do%20files/Scenario%20A%20-%20logit%20-%20runs%20sim%20study%20I%20for%20Monte%20Carlo%20QBA.do):
            Runs scenario A of simulation study I for Monte Carlo QBA.
        - [**Results/**](./Simulation%20study/Simulation%20study%20I/Scenario%20A/Monte%20Carlo/Results/)
          - empty folder to store the results from running the above do file
      - [**Bayesian/**](./Simulation%20study/Simulation%20study%20I/Scenario%20A/Bayesian/)
        - [**Seeds/**](./Simulation%20study/Simulation%20study%20I/Scenario%20A/Bayesian/Seeds/)
          - [Generates random
            integers.do](./Simulation%20study/Simulation%20study%20I/Scenario%20A/Bayesian/Seeds/Generates random integers.do):
            Generates 8,000 random integers for seeds.
          - [RandomIntegers.csv](./Simulation%20study/Simulation%20study%20I/Scenario%20A/Bayesian/Seeds/RandomIntegers.csv):
            Lists 8,000 random integers.
        - [**R
          files/**](./Simulation%20study/Simulation%20study%20I/Scenario%20A/Bayesian/R%20files/)
          - [ScenarioA-glm-RunsSimStudyIforBayesianQBA.R](./Simulation%20study/Simulation%20study%20I/Scenario%20A/Bayesian/R%20files/ScenarioA-glm-RunsSimStudyIforBayesianQBA.R) :
            Runs scenario A of simulation study I for Bayesian QBA.
        - [**Results/**](./Simulation%20study/Simulation%20study%20I/Scenario%20A/Bayesian/Results/)
          - empty folder to store the results from running the above R file
    - [**Scenario
      B/**](./Simulation%20study/Simulation%20study%20I/Scenario%20B/)
      - [**Data/**](./Simulation%20study/Simulation%20study%20I/Scenario%20B/Data/)
        - 500 csv files from
          [Dataset_1.csv](./Simulation%20study/Simulation%20study%20I/Scenario%20B/Data/Dataset_1.csv)
          to
          [Dataset_500.csv](./Simulation%20study/Simulation%20study%20I/Scenario%20B/Data/Dataset_500.csv)
          of the simulated datasets for scenario B of simulation studies I, II
          and VI
      - [**Monte
        Carlo/**](./Simulation%20study/Simulation%20study%20I/Scenario%20B/Monte%20Carlo)
        - [**Do
          files/**](./Simulation%20study/Simulation%20study%20I/Scenario%20B/Monte%20Carlo/Do%20files/)
          - [Scenario B - regress - runs sim study I for Monte Carlo
            QBA.do](./Simulation%20study/Simulation%20study%20I/Scenario%20B/Monte%20Carlo/Do%20files/Scenario%20B%20-%20regress%20-%20runs%20sim%20study%20I%20for%20Monte%20Carlo%20QBA.do):
            Runs scenario B of simulation study I for Monte Carlo QBA.
        - [**Results/**](./Simulation%20study/Simulation%20study%20I/Scenario%20B/Monte%20Carlo/Results/)
          - empty folder to store the results from running the above do file
      - [**Bayesian/**](./Simulation%20study/Simulation%20study%20I/Scenario%20B/Bayesian)
        - [**Seeds/**](./Simulation%20study/Simulation%20study%20I/Scenario%20B/Bayesian/Seeds/)
          - [Generates random
            integers.do](./Simulation%20study/Simulation%20study%20I/Scenario%20B/Bayesian/Seeds/Generates random integers.do):
            Generates 8,000 random integers for seeds.
          - [RandomIntegers.csv](./Simulation%20study/Simulation%20study%20I/Scenario%20B/Bayesian/Seeds/RandomIntegers.csv):
            Lists 8,000 random integers.
        - [**R
          files/**](./Simulation%20study/Simulation%20study%20I/Scenario%20B/Bayesian/R%20files/)
          - [ScenarioB-lm-RunsSimStudyIforBayesianQBA.R](./Simulation%20study/Simulation%20study%20I/Scenario%20B/Bayesian/R%20files/ScenarioB-lm-RunsSimStudyIforBayesianQBA.R):
            Runs scenario B of simulation study I for Bayesian QBA.
        - [**Results/**](./Simulation%20study/Simulation%20study%20I/Scenario%20B/Bayesian/Results/)
          - empty folder to store the results from running the above R file

  - [**Simulation study II/**](./Simulation%20study/Simulation%20study%20II/)
    - [**Do files/**](./Simulation%20study/Simulation%20study%20II/Do%20files/)
      - [Scenario A - logit - runs sim study II for Monte Carlo
        QBA.do](./Simulation%20study/Simulation%20study%20II/Do%20files/Scenario A - logit - runs sim study II for Monte Carlo QBA.do):
        Runs scenario A of simulation study II for Monte Carlo QBA.
      - [Scenario B - regress - runs sim study II for Monte Carlo
        QBA.do](./Simulation%20study/Simulation%20study%20II/Do%20files/Scenario B - regress - runs sim study II for Monte Carlo QBA.do):
        Runs scenario B of simulation study II for Monte Carlo QBA.
    - [**Results/**](./Simulation%20study/Simulation%20study%20II/Results/)
      - empty folder to store the results from running the above do files

  - [**Simulation study III/**](./Simulation%20study/Simulation%20study%20III/)
    - [**Do files/**](./Simulation%20study/Simulation%20study%20III/Do%20files/)
      - [Scenario C - logit - runs sim study III for Monte Carlo
        QBA.do](./Simulation%20study/Simulation%20study%20III/Do%20files/Scenario C - logit - runs sim study III for Monte Carlo QBA.do):
        Runs scenario C of simulation study III for Monte Carlo QBA.
      - [Scenario D - regress - runs sim study III for Monte Carlo
        QBA.do](./Simulation%20study/Simulation%20study%20III/Do%20files/Scenario D - regress - runs sim study III for Monte Carlo QBA.do):
        Runs scenario D of simulation study III for Monte Carlo QBA.
    - [**Results/**](./Simulation%20study/Simulation%20study%20III/Results/)
      - empty folder to store the results from running the above do files

  - [**Simulation study IV/**](./Simulation%20study/Simulation%20study%20IV/)
    - [**Data/**](./Simulation%20study/Simulation%20study%20IV/Data/)
      - 500 csv files from
        [Dataset_1.csv](./Simulation%20study/Simulation%20study%20IV/Data/Dataset_1.csv)
        to
        [Dataset_500.csv](./Simulation%20study/Simulation%20study%20IV/Data/Dataset_500.csv)
        of the simulated datasets for scenario G of simulation study IV
    - [**Do files/**](./Simulation%20study/Simulation%20study%20IV/Do%20files/)
      - [Scenario E - regress - runs sim study IV for Monte Carlo
        QBA.do](./Simulation%20study/Simulation%20study%20IV/Do%20files/Scenario%20E%20-%20regress%20-%20runs%20sim%20study%20IV%20for%20Monte%20Carlo%20QBA.do):
        Runs scenario E of simulation study IV for Monte Carlo QBA.
      - [Scenario F - mlogit - runs sim study IV for Monte Carlo
        QBA.do](./Simulation%20study/Simulation%20study%20IV/Do%20files/Scenario F - mlogit - runs sim study IV for Monte Carlo QBA.do):
        Runs scenario F of simulation study IV for Monte Carlo QBA.
      - [Scenario G - Cox PH - runs sim study IV for Monte Carlo
        QBA.do](./Simulation%20study/Simulation%20study%20IV/Do%20files/Scenario%20G%20-%20Cox%20PH%20-%20runs%20sim%20study%20IV%20for%20Monte%20Carlo%20QBA.do):
        Runs scenario G of simulation study IV for Monte Carlo QBA.
    - [**Results/**](./Simulation%20study/Simulation%20study%20IV/Results/)
      - empty folder to store the results from running the above do files

  - [**Simulation study V/**](./Simulation%20study/Simulation%20study%20V/)
    - [**Data**](./Simulation%20study/Simulation%20study%20V/Data/)
      - 500 dta files from
        [Dataset_1.dta](./Simulation%20study/Simulation%20study%20V/Data/Dataset_1.dta)
        to
        [Dataset_500.dta](./Simulation%20study/Simulation%20study%20V/Data/Dataset_500.dta)
        of the simulated datasets for scenarios A and B of simulation study V
    - [**Do files/**](./Simulation%20study/Simulation%20study%20V/Do%20files/)
      - [Scenario A - logit - runs sim study V for Monte Carlo
        QBA.do](./Simulation%20study/Simulation%20study%20V/Do%20files/Scenario%20A%20-%20logit%20-%20runs%20sim%20study%20V%20for%20Monte%20Carlo%20QBA.do):
        Runs scenario A of simulation study V for Monte Carlo QBA.
      - [Scenario B - regress - runs sim study V for Monte Carlo
        QBA.do](./Simulation%20study/Simulation%20study%20V/Do%20files/Scenario%20B%20-%20regress%20-%20runs%20sim%20study%20V%20for%20Monte%20Carlo%20QBA.do):
        Runs scenario B of simulation study V for Monte Carlo QBA.
    - [**Results/**](./Simulation%20study/Simulation%20study%20V/Results/)
      - empty folder to store the results from running the above do files

  - [**Simulation study VI/**](./Simulation%20study/Simulation%20study%20VI/)
    - [**Do files/**](./Simulation%20study/Simulation%20study%20VI/Do%20files/)
      - [Scenario A - logit - runs sim study VI for Monte Carlo
        QBA.do](./Simulation%20study/Simulation%20study%20VI/Do%20files/Scenario%20A%20-%20logit%20-%20runs%20sim%20study%20VI%20for%20Monte%20Carlo%20QBA.do):
        Runs scenario A of simulation study VI for Monte Carlo QBA.
      - [Scenario B - regress - runs sim study VI for Monte Carlo
        QBA.do](./Simulation%20study/Simulation%20study%20VI/Do%20files/Scenario%20B%20-%20regress%20-%20runs%20sim%20study%20VI%20for%20Monte%20Carlo%20QBA.do):
        Runs scenario B of simulation study VI for Monte Carlo QBA.
    - [**Results/**](./Simulation%20study/Simulation%20study%20VI/Results/)
      - empty folder to store the results from running the above do files

- [**Applied example/**](./Applied%20example/)
  - [GeneratesDataForNHANESExample.R](./Applied%20example/GeneratesDataForNHANESExample.R):
    Generates data for the NHANES example.
  - [**Data/**](./Applied%20example/Data/)
    - [Data.csv](./Applied%20example/Data/Data.csv): cleaned NHANES data for
      analysis
  - [**Monte Carlo/**](./Applied%20example/Monte%20Carlo/)
    - [Apply MCQBA to NHANES example using prior for
      piM.do](./Applied%20example/Monte%20Carlo/Apply%20MCQBA%20to%20NHANES %20%20%20%20%20%20example%20using%20prior%20for%20piM.do):
      Applies Monte Carlo QBA with a prior distribution for the marginal
      prevalence.
    - [Apply MCQBA to NHANES example using prior for the
      intercept.do](./Applied%20example/Monte%20Carlo/Apply%20MCQBA%20to%20NHANES%20example%20using%20prior%20for%20the%20intercept.do):
      Applies Monte Carlo QBA with a prior distribution for the intercept of
      logistic regression.
    - [Apply MCQBA to NHANES example using inaccurate
      priors.do](./Applied%20example/Monte%20Carlo/Apply%20MCQBA%20to%20NHANES%20example%20using%20inaccurate%20priors.do):
      Applies Monte Carlo QBA using inaccurate priors.
  - [**Bayesian/**](./Applied%20example/Bayesian/)
    - [ApplyBayesianQBAToNHANESExample.R](./Applied%20example/Bayesian/ApplyBayesianQBAToNHANESExample.R):
      Applies Bayesian QBA.

- [**Bayesian function/**](./Bayesian%20function/)
  - [MyUnm_glm.R](./Bayesian%20function/MyUnm_glm.R): modified *unm_glm*
    function from R package *unmconf* 1.0.0 for running Bayesian QBA

## Publications

Emily Kawabata, Chin Yang Shapland, Tom Palmer, David Carslake, Kate Tilling,
Rachael Hughes. A flexible Monte Carlo quantitative bias analysis for unmeasured
confounding. Statistical Methods in Medical Research (accepted in 2026).

See also:
[https://www.medrxiv.org/content/10.1101/2025.08.12.25333217v1](https://www.medrxiv.org/content/10.1101/2025.08.12.25333217v1)

## Further information

If you would like any further information, please contact
emily.kawabata@bristol.ac.uk.
