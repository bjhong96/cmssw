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

echo "#!/bin/sh" 							   >  ${dirname}/submit_all.sh
echo " " 							  	   >> ${dirname}/submit_all.sh
echo "for i in {1..${linenum_tot}}; do" 				   >> ${dirname}/submit_all.sh
echo "  config_file=\"muonValidation_${p_or_np}_${linenum_tot}_\${i}.sub\"" >> ${dirname}/submit_all.sh
echo "  condor_submit \${config_file}" 					   >> ${dirname}/submit_all.sh
echo "done" 								   >> ${dirname}/submit_all.sh

chmod +x ${dirname}/submit_all.sh
