# Pycnandra_V2
Pycnandra_V2

### Installing pipeline :

First, open your terminal. Then, run these two command lines :

    cd -place_in_your_local_computer
    git clone https://github.com/PLStenger/Pycnandra_V2.git

### Run this :

    time nohup bash 001_pipeline_QIIME2_PE_Pycnandra_V2.sh &>  001_pipeline_QIIME2_PE_Pycnandra_V2.out
    time nohup bash 002_rarefaction_Pycnandra_V2.sh &> 002_rarefaction_Pycnandra_V2.out
    time nohup bash 003_qiime2_assign_taxonomy_Pycnandra_V2.sh &> 003_qiime2_assign_taxonomy_Pycnandra_V2.out
    time nohup bash 004_core_biom_Pycnandra_V2.sh &> 004_core_biom_Pycnandra_V2.out
