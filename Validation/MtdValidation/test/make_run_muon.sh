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

dir_src=$(echo "$PWD" | grep -o '.*/src')

while read line; do
  ((linenum+=1))
  echo $linenum

  echo "#!/bin/bash" 				>  ${dirname}/run_${p_or_np}_${linenum_tot}_${linenum}.sh
  echo " " 					>> ${dirname}/run_${p_or_np}_${linenum_tot}_${linenum}.sh
  echo "echo \"start\"" 			>> ${dirname}/run_${p_or_np}_${linenum_tot}_${linenum}.sh
  echo "echo \"\$PWD\"" 			>> ${dirname}/run_${p_or_np}_${linenum_tot}_${linenum}.sh
  echo "cd $dir_src" 				>> ${dirname}/run_${p_or_np}_${linenum_tot}_${linenum}.sh
  echo "cmsenv" 				>> ${dirname}/run_${p_or_np}_${linenum_tot}_${linenum}.sh
  echo "sleep 30s" 				>> ${dirname}/run_${p_or_np}_${linenum_tot}_${linenum}.sh
  echo "cd $PWD/$dirname" 			>> ${dirname}/run_${p_or_np}_${linenum_tot}_${linenum}.sh
  echo "echo \"\$PWD\"" 			>> ${dirname}/run_${p_or_np}_${linenum_tot}_${linenum}.sh
  echo "cmsRun mtdValidation_muon_cfg_${p_or_np}_${linenum_tot}_${linenum}.py" >> ${dirname}/run_${p_or_np}_${linenum_tot}_${linenum}.sh

  done < ${dirname}/flist.txt
