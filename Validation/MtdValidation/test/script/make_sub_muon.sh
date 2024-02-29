#!/bin/sh

dirname=$1
#dirname="240227_PU200_prompt_muon_test"
mkdir -p ${dirname}
mkdir -p ${dirname}/output
mkdir -p ${dirname}/log

if [[ "$dirname" == *"_prompt"* ]]; then
  echo "samples for prompt muon"
  p_or_np="prompt"
elif [[ "$dirname" == *"_nonprompt"* ]]; then
  echo "samples for nonprompt muon"
  p_or_np="nonprompt"
else
  echo "dirname should include \"prompt\" or \"nonprompt\""
  break
fi

if [[ "$dirname" == *"_PU200"* ]]; then
  echo "samples for PU200 scenario"
  pu_or_nopu="PU200"
elif [[ "$dirname" == *"_noPU"* ]]; then
  echo "samples for noPU scenario"
  pu_or_nopu="noPU"
else
  echo "dirname should include \"PU200\" or \"noPU\""
  break
fi

linenum_tot=$(awk "END { print NR }" ${dirname}/flist.txt)
linenum=0

while read line; do
  ((linenum+=1))
  echo $linenum

  echo "batch_name = muonValidation_${pu_or_nopu}_${p_or_np}_${linenum_tot}_${linenum}"  			  > ${dirname}/muonValidation_${p_or_np}_${linenum_tot}_${linenum}.sub
  echo "executable = run_${p_or_np}_${linenum_tot}_${linenum}.sh" 						  >> ${dirname}/muonValidation_${p_or_np}_${linenum_tot}_${linenum}.sub
  echo "universe   = vanilla" 											  >> ${dirname}/muonValidation_${p_or_np}_${linenum_tot}_${linenum}.sub
  echo "output                          = log/muonValidation_${p_or_np}_${linenum_tot}_${linenum}_\$(ProcId).out" >> ${dirname}/muonValidation_${p_or_np}_${linenum_tot}_${linenum}.sub
  echo "error                           = log/muonValidation_${p_or_np}_${linenum_tot}_${linenum}_\$(ProcId).err" >> ${dirname}/muonValidation_${p_or_np}_${linenum_tot}_${linenum}.sub
  echo "log                             = log/muonValidation_${p_or_np}_${linenum_tot}_${linenum}_\$(ProcId).log" >> ${dirname}/muonValidation_${p_or_np}_${linenum_tot}_${linenum}.sub
  echo "should_transfer_files           = YES"           >> ${dirname}/muonValidation_${p_or_np}_${linenum_tot}_${linenum}.sub
  echo "when_to_transfer_output         = ON_EXIT"       >> ${dirname}/muonValidation_${p_or_np}_${linenum_tot}_${linenum}.sub
  echo "+JobFlavour                     = \"testmatch\"" >> ${dirname}/muonValidation_${p_or_np}_${linenum_tot}_${linenum}.sub
  echo " " 						 >> ${dirname}/muonValidation_${p_or_np}_${linenum_tot}_${linenum}.sub
  echo "RequestMemory = 20000" 				 >> ${dirname}/muonValidation_${p_or_np}_${linenum_tot}_${linenum}.sub
  echo "RequestCpus   = 4" 				 >> ${dirname}/muonValidation_${p_or_np}_${linenum_tot}_${linenum}.sub
  echo " " 						 >> ${dirname}/muonValidation_${p_or_np}_${linenum_tot}_${linenum}.sub
  echo "Queue" 						 >> ${dirname}/muonValidation_${p_or_np}_${linenum_tot}_${linenum}.sub

done < ${dirname}/flist.txt
