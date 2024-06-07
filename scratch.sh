# run all linear regression

# 10sp
./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_10sp_random_1000 1000 0 100 0 1
./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_10sp_random_1000 1000 0 10000 0 1

./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_10sp_random_1000 1000 1 100 0 1
./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_10sp_random_1000 1000 1 10000 0 1

# for test set
./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_10sp_random_testSet 100 0 100 0 1
./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_10sp_random_testSet 100 0 10000 0 1

./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_10sp_random_testSet 100 1 100 0 1
./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_10sp_random_testSet 100 1 10000 0 1

# 10sp
./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_10sp 100 0 100 0 1
./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_10sp 100 0 10000 0 1

./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_10sp 100 1 100 0 1
./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_10sp 100 1 10000 0 1

#logi
./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_10sp 100 0 100 1 0
./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_10sp 100 0 10000 1 0

./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_10sp 100 1 100 1 0
./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_10sp 100 1 10000 1 0

# 100sp
./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_100sp 100 0 100 0 1
./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_100sp 100 0 1000 0 1
./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_100sp 100 0 10000 0 1

./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_100sp 100 1 100 0 1
./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_100sp 100 1 1000 0 1
./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_100sp 100 1 10000 0 1

#logi
./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_100sp 100 0 100 1 0
./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_100sp 100 0 10000 1 0

./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_100sp 100 1 100 1 0
./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_100sp 100 1 10000 1 0


# 10sp_random
./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_10sp_random_moreSamples 100 0 100 0 1
./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_10sp_random_moreSamples 100 0 10000 0 1

./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_10sp_random_moreSamples 100 1 100 0 1
./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_10sp_random_moreSamples 100 1 10000 0 1

# 100sp_random
./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_100sp_random_moreSamples 100 0 100 0 1
./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_100sp_random_moreSamples 100 0 10000 0 1

./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_100sp_random_moreSamples 100 1 100 0 1
./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/multiSim_100sp_random_moreSamples 100 1 10000 0 1
