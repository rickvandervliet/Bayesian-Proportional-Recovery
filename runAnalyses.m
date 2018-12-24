%% Step 1. Proportional Recovery Overfit Using K=10
disp('%% Step 1. Proportional Recovery Overfit Using K=10')
proportionalRecoveryOverfit

%% Step 2. Refit data using Ktrue
disp('%% Step 2. Refit data using Ktrue')
proportionalRecoveryFit

%% Step 3. Cross-validation for prediction
disp('%% Step 3. Cross-validation for prediction')
proportionalRecoveryPrediction

%% Step 4. Plot results
disp('%% Step 4. Plot results')
resultsForPaper