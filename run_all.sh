#!/bin/bash
LOG="$HOME/full_scale_experiments.log"

STEP1="addpath('experiments');run_comparison_para('n_runs',30,'map_name','Map1_Small')"
STEP2="addpath('experiments');run_comparison_para('n_runs',30,'map_name','Map1_Medium')"
STEP3="addpath('experiments');run_comparison_para('n_runs',30,'map_name','Map1_Large')"
STEP4="addpath('experiments');run_ablation_para('n_runs',30,'map_name','Map1_Small')"
STEP5="addpath('experiments');run_ablation_para('n_runs',30,'map_name','Map1_Medium')"
STEP6="addpath('experiments');run_ablation_para('n_runs',30,'map_name','Map1_Large')"

for STEP in "$STEP1" "$STEP2" "$STEP3" "$STEP4" "$STEP5" "$STEP6"
do
  echo "========== Running: $STEP ==========" | tee -a $LOG
  cd ~/csa_goa3clear1
  matlab -batch "cp=parcluster('local');cp.NumWorkers=4;parpool(cp,4);$STEP;delete(gcp('nocreate'));" 2>&1 | tee -a $LOG
  if [ $? -ne 0 ]; then
    echo "WARNING: 失败，继续执行后续实验..." | tee -a $LOG
  fi
  echo "" | tee -a $LOG
done
echo "========== 全部完成 ==========" | tee -a $LOG
