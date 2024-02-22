#Fisher exact tests for table 1 and table 2 entries where the Chi-squared test wasn't appropriate due to less than 5 observation in a category



### TABLE 2 ###


## THC freq ##

#contingency table
THC_freq_contingency_table <- matrix(c(27, 6, 2, 2, 6, 4, 2, 9, 4, 9, 0, 3, 0, 8), nrow = 7, byrow = TRUE)

# Perform Fisher's exact test
THC_freq_fisher_result <- fisher.test(THC_freq_contingency_table)

# Print results
print(THC_freq_fisher_result)




### TABLE 1 + paired ###

## race + paired ##

#contingency table
race_paired_contingency_table <- matrix(c(3,3,1,1,2,6,4,2,27,45,48,37,0,2,1,1,0,1,0,0), nrow = 5, byrow = TRUE)

# Perform Fisher's exact test
race_paired_fisher_result <- fisher.test(race_paired_contingency_table)

# Print results
print(race_paired_fisher_result)



## hispanic + paired ##

#contingency table
hispanic_paired_contingency_table <- matrix(c(3,1,1,1,29,56,53,40), nrow = 2, byrow = TRUE)

# Perform Fisher's exact test
hispanic_paired_fisher_result <- fisher.test(hispanic_paired_contingency_table)

# Print results
print(hispanic_paired_fisher_result)



## education + paired ##

#contingency table
education_paired_contingency_table <- matrix(c(0,3,2,1,2,9,5,5,1,2,1,1,10,18,18,12,19,25,28,22), nrow = 5, byrow = TRUE)

# Perform Fisher's exact test
education_paired_fisher_result <- fisher.test(education_paired_contingency_table,workspace=2e7)

# Print results
print(education_paired_fisher_result)



## employment + paired ##

#contingency table
employment_paired_contingency_table <- matrix(c(24,29,31,22,2,5,5,4,1,1,0,0,5,11,8,6,0,4,2,2,0,1,1,1,0,2,3,2,0,4,4,4), nrow = 8, byrow = TRUE)

# Perform Fisher's exact test
employment_paired_fisher_result <- fisher.test(employment_paired_contingency_table,workspace=2e7)

# Print results
print(employment_paired_fisher_result)



## handedness + paired ##

#contingency table
handedness_paired_contingency_table <- matrix(c(3,8,8,6,29,49,46,35), nrow = 2, byrow = TRUE)

# Perform Fisher's exact test
handedness_paired_fisher_result <- fisher.test(handedness_paired_contingency_table)

# Print results
print(handedness_paired_fisher_result)




### TABLE 1 + unpaired ###

## race + unpaired ##

#contingency table
race_unpaired_contingency_table <- matrix(c(3,3,1,2,6,4,27,45,48,0,2,1,0,1,0), nrow = 5, byrow = TRUE)

# Perform Fisher's exact test
race_unpaired_fisher_result <- fisher.test(race_unpaired_contingency_table)

# Print results
print(race_unpaired_fisher_result)



## hispanic + unpaired ##

#contingency table
hispanic_unpaired_contingency_table <- matrix(c(3,1,1,29,56,53), nrow = 2, byrow = TRUE)

# Perform Fisher's exact test
hispanic_unpaired_fisher_result <- fisher.test(hispanic_unpaired_contingency_table)

# Print results
print(hispanic_unpaired_fisher_result)



## education + unpaired ##

#contingency table
education_unpaired_contingency_table <- matrix(c(0,3,2,2,9,5,1,2,1,10,18,18,19,25,28), nrow = 5, byrow = TRUE)

# Perform Fisher's exact test
education_unpaired_fisher_result <- fisher.test(education_unpaired_contingency_table,workspace=2e7)

# Print results
print(education_unpaired_fisher_result)



## employment + unpaired ##

#contingency table
employment_unpaired_contingency_table <- matrix(c(24,29,31,2,5,5,1,1,0,5,11,8,0,4,2,0,1,1,0,2,3,0,4,4), nrow = 8, byrow = TRUE)

# Perform Fisher's exact test
employment_unpaired_fisher_result <- fisher.test(employment_unpaired_contingency_table,workspace=2e7)

# Print results
print(employment_unpaired_fisher_result)



## handedness + unpaired ##

#contingency table
handedness_unpaired_contingency_table <- matrix(c(3,8,8,29,49,46), nrow = 2, byrow = TRUE)

# Perform Fisher's exact test
handedness_unpaired_fisher_result <- fisher.test(handedness_unpaired_contingency_table)

# Print results
print(handedness_unpaired_fisher_result)
