# Bayesian-Proportional-Recovery

The files in this GitHub accompany the papers:
1. Predicting recovery of the upper limb after stroke; a longitudinal mixture model
2. How to improve the study power of upper limb recovery and rehabilitation trials early post stroke?

################################################################################################################
1. Predicting recovery of the upper limb after stroke; a longitudinal mixture model
In the first paper, we developed an longitudinal exponential recovery mixture model for upper extremity recovery 
on the Fugl-Meyer upper extremity assessment in the subacute phase after stroke.
The following files were used for the analysis in this paper:

1. proportionalRecoveryOverfit.m, and proportionalRecovery.txt
These files implement the longitudinal mixture model. The number of groups is determined based on the Rousseau 
and Mengerson criterion by choosing a large number of subgroups (10) and selecting only those subgroups which
have a certain assignment probability.

2. clusterBugs.m, proportionalRecoveryCluster.m, proportionalRecoveryCluster.txt
These files are used to assign a stroke recovery cluster to all patients who were initially excluded from the 
analysis because they deteriorated clinically. This stroke recovery cluster is used later in the 
cross-validation step.

3. proportionalRecoveryPrediction.m, fitBugs.m, proportionalRecovery.txt, predictionBugs.m, 
proportionalRecoveryPrediction.txt
These files are used for cross-validation of the predictions. First, a new model is created using data from all
but n-1 patients. This model is then applied to sequential data of the remaining patient to predict outcome.

################################################################################################################
2. How to improve the study power of upper limb recovery and rehabilitation trials early post stroke?
In the second paper, we amended the model of the first paper with an extra term for an intervention effect.

