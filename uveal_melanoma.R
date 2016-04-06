#untar('gdac.broadinstitute.org_UVM.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz', exdir = '.')

#untar('gdac.broadinstitute.org_UVM.Merge_Clinical.Level_1.2016012800.0.0.tar.gz', exdir = '.')

clinical <-  read.delim('./gdac.broadinstitute.org_UVM.Merge_Clinical.Level_1.2016012800.0.0/UVM.merged_only_clinical_clin_format.txt',
                        skip = 1)
clinical <-  read.delim('./gdac.broadinstitute.org_UVM.Clinical_Pick_Tier1.Level_4.2016012800.0.0/UVM.clin.merged.picked.txt')


