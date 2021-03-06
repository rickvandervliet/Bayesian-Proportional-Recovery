# Bayesian-Proportional-Recovery

The files in this GitHub accompany the papers:
1. Predicting recovery from upper limb motor impairment after stroke: a longitudinal mixture model
2. Improving power of upper limb motor rehabilitation trials with a longitudinal mixture model

## Predicting recovery from upper limb motor impairment after stroke: a longitudinal mixture model

In the first paper, we developed an longitudinal exponential recovery mixture model for upper extremity recovery 
on the Fugl-Meyer upper extremity assessment in the subacute phase after stroke.
The following files were used for the analysis in this paper:

1. proportionalRecoveryOverfit.m, and proportionalRecovery.txt
These files implement the longitudinal mixture model. The number of groups is determined based on the Rousseau 
and Mengerson criterion by choosing a large number of subgroups (10) and selecting only those subgroups which
have a certain assignment probability.

2. proportionalRecoveryPrediction.m, fitBugs.m, proportionalRecovery.txt, predictionBugs.m, 
proportionalRecoveryPrediction.txt
These files are used for cross-validation of the predictions. First, a new model is created using data from all
but n-1 patients. This model is then applied to sequential data of the remaining patient to predict outcome.

Clinicians can make use of the following matlab files to make predictions for their patients:
- clinicalPredicationsFMUE.m, which predicts the FM-UE at weekly time points up until 26 weeks post stroke
- clinicalPredictionsFMUEplot.m, which plots the time course and the histogram of the FM-UE at 26 weeks post stroke

## Improving power of upper limb motor rehabilitation trials with a longitudinal mixture model

In the second paper, we amended the model of the first paper with an extra term for an intervention effect.
This way, the model can be used to analyse stroke recovery and rehabilitation trials.

The following files were used for the analysis in this paper:

1. generateBugs.m and proportionalRecoveryGeneration.txt
To create simulated datasets based on the model parameters derived in the first paper.

2. powerBugs.m and proportionalRecoveryPower.txt 
To estimate the intervention effect from a simulated study

3. proportionalRecoveryPower.m
Master script to generate data and estimate intervention effects for different study designs.