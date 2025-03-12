# 26Mg_an_analysis
 26Mg(a,n) analysis
 To convert binary file to root file: ./scan run_1525
 
 To compile the unfolding code: g++ Unfolding_MLEM_covar.cpp -I. -o Unfolding_MLEM_25Mg.exe
To run the unfolding code: root -l Analysis_singleChn_25Mg_new.cpp
root -l Analysis_singleChn_25Mg_new_neutron_energies.cpp
root -l Analysis_singleChn.cpp
