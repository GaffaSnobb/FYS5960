#!/bin/sh 
# export OMP_STACKSIZE=1g
export GFORTRAN_UNBUFFERED_PRECONNECTED=y
# ulimit -s unlimited

# ---------- Ni68_jun45 --------------
echo "start running log_Ni68_jun45_m0p.txt ..."
cat > Ni68_jun45_0.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 1.1
  fn_int = "jun45.snt"
  fn_ptn = "Ni68_jun45_p.ptn"
  fn_save_wave = "Ni68_jun45_m0p.wav"
  gl = 1.0, 0.0
  gs = 3.910, -2.678
  hw_type = 1
  is_double_j = .false.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 0
  n_block = 0
  n_eigen = 10
  n_restart_vec = 15
  orbs_ratio = 4, 8
&end
EOF
nice ./kshell.exe Ni68_jun45_0.input > log_Ni68_jun45_m0p.txt 2>&1 

rm -f tmp_snapshot_Ni68_jun45_p.ptn_0_* tmp_lv_Ni68_jun45_p.ptn_0_* Ni68_jun45_0.input 


# --------------- transition probabilities --------------

echo "start running log_Ni68_jun45_tr_m0p_m0p.txt ..."
cat > Ni68_jun45_0_0.input <<EOF
&input
  fn_int   = "jun45.snt"
  fn_ptn_l = "Ni68_jun45_p.ptn"
  fn_ptn_r = "Ni68_jun45_p.ptn"
  fn_load_wave_l = "Ni68_jun45_m0p.wav"
  fn_load_wave_r = "Ni68_jun45_m0p.wav"
  hw_type = 1
  eff_charge = 1.5, 1.1
  gl = 1.0, 0.0
  gs = 3.910, -2.678
&end
EOF
nice ./transit.exe Ni68_jun45_0_0.input > log_Ni68_jun45_tr_m0p_m0p.txt 2>&1 

rm -f Ni68_jun45_0_0.input


./collect_logs.py log_*Ni68_jun45* > summary_Ni68_jun45.txt
echo "Finish computing Ni68_jun45.    See summary_Ni68_jun45.txt"
echo 

