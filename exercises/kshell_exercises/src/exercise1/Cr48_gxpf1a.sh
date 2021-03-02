#!/bin/sh 
# export OMP_STACKSIZE=1g
export GFORTRAN_UNBUFFERED_PRECONNECTED=y
# ulimit -s unlimited

# ---------- Cr48_gxpf1a --------------
echo "start running log_Cr48_gxpf1a_j0p.txt ..."
cat > Cr48_gxpf1a_0.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5
  fn_int = "gxpf1a.snt"
  fn_ptn = "Cr48_gxpf1a_p.ptn"
  fn_save_wave = "Cr48_gxpf1a_j0p.wav"
  gl = 1.0, 0.0, 
  gs = 3.91, -2.678, 
  hw_type = 1
  is_double_j = .true.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 0
  n_block = 0
  n_eigen = 1
  n_restart_vec = 10
  orbs_ratio = 2, 3, 4, 6, 7, 8
&end
EOF
nice ./kshell.exe Cr48_gxpf1a_0.input > log_Cr48_gxpf1a_j0p.txt 2>&1 

rm -f tmp_snapshot_Cr48_gxpf1a_p.ptn_0_* tmp_lv_Cr48_gxpf1a_p.ptn_0_* Cr48_gxpf1a_0.input 


echo "start running log_Cr48_gxpf1a_j4p.txt ..."
cat > Cr48_gxpf1a_4.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5
  fn_int = "gxpf1a.snt"
  fn_ptn = "Cr48_gxpf1a_p.ptn"
  fn_save_wave = "Cr48_gxpf1a_j4p.wav"
  gl = 1.0, 0.0, 
  gs = 3.91, -2.678, 
  hw_type = 1
  is_double_j = .true.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 4
  n_block = 0
  n_eigen = 1
  n_restart_vec = 10
  orbs_ratio = 2, 3, 4, 6, 7, 8
&end
EOF
nice ./kshell.exe Cr48_gxpf1a_4.input > log_Cr48_gxpf1a_j4p.txt 2>&1 

rm -f tmp_snapshot_Cr48_gxpf1a_p.ptn_4_* tmp_lv_Cr48_gxpf1a_p.ptn_4_* Cr48_gxpf1a_4.input 


echo "start running log_Cr48_gxpf1a_j8p.txt ..."
cat > Cr48_gxpf1a_8.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5
  fn_int = "gxpf1a.snt"
  fn_ptn = "Cr48_gxpf1a_p.ptn"
  fn_save_wave = "Cr48_gxpf1a_j8p.wav"
  gl = 1.0, 0.0, 
  gs = 3.91, -2.678, 
  hw_type = 1
  is_double_j = .true.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 8
  n_block = 0
  n_eigen = 1
  n_restart_vec = 10
  orbs_ratio = 2, 3, 4, 6, 7, 8
&end
EOF
nice ./kshell.exe Cr48_gxpf1a_8.input > log_Cr48_gxpf1a_j8p.txt 2>&1 

rm -f tmp_snapshot_Cr48_gxpf1a_p.ptn_8_* tmp_lv_Cr48_gxpf1a_p.ptn_8_* Cr48_gxpf1a_8.input 


echo "start running log_Cr48_gxpf1a_j12p.txt ..."
cat > Cr48_gxpf1a_12.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5
  fn_int = "gxpf1a.snt"
  fn_ptn = "Cr48_gxpf1a_p.ptn"
  fn_save_wave = "Cr48_gxpf1a_j12p.wav"
  gl = 1.0, 0.0, 
  gs = 3.91, -2.678, 
  hw_type = 1
  is_double_j = .true.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 12
  n_block = 0
  n_eigen = 1
  n_restart_vec = 10
  orbs_ratio = 2, 3, 4, 6, 7, 8
&end
EOF
nice ./kshell.exe Cr48_gxpf1a_12.input > log_Cr48_gxpf1a_j12p.txt 2>&1 

rm -f tmp_snapshot_Cr48_gxpf1a_p.ptn_12_* tmp_lv_Cr48_gxpf1a_p.ptn_12_* Cr48_gxpf1a_12.input 


echo "start running log_Cr48_gxpf1a_j16p.txt ..."
cat > Cr48_gxpf1a_16.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5
  fn_int = "gxpf1a.snt"
  fn_ptn = "Cr48_gxpf1a_p.ptn"
  fn_save_wave = "Cr48_gxpf1a_j16p.wav"
  gl = 1.0, 0.0, 
  gs = 3.91, -2.678, 
  hw_type = 1
  is_double_j = .true.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 16
  n_block = 0
  n_eigen = 1
  n_restart_vec = 10
  orbs_ratio = 2, 3, 4, 6, 7, 8
&end
EOF
nice ./kshell.exe Cr48_gxpf1a_16.input > log_Cr48_gxpf1a_j16p.txt 2>&1 

rm -f tmp_snapshot_Cr48_gxpf1a_p.ptn_16_* tmp_lv_Cr48_gxpf1a_p.ptn_16_* Cr48_gxpf1a_16.input 


echo "start running log_Cr48_gxpf1a_j20p.txt ..."
cat > Cr48_gxpf1a_20.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5
  fn_int = "gxpf1a.snt"
  fn_ptn = "Cr48_gxpf1a_p.ptn"
  fn_save_wave = "Cr48_gxpf1a_j20p.wav"
  gl = 1.0, 0.0, 
  gs = 3.91, -2.678, 
  hw_type = 1
  is_double_j = .true.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 20
  n_block = 0
  n_eigen = 1
  n_restart_vec = 10
  orbs_ratio = 2, 3, 4, 6, 7, 8
&end
EOF
nice ./kshell.exe Cr48_gxpf1a_20.input > log_Cr48_gxpf1a_j20p.txt 2>&1 

rm -f tmp_snapshot_Cr48_gxpf1a_p.ptn_20_* tmp_lv_Cr48_gxpf1a_p.ptn_20_* Cr48_gxpf1a_20.input 


echo "start running log_Cr48_gxpf1a_j24p.txt ..."
cat > Cr48_gxpf1a_24.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5
  fn_int = "gxpf1a.snt"
  fn_ptn = "Cr48_gxpf1a_p.ptn"
  fn_save_wave = "Cr48_gxpf1a_j24p.wav"
  gl = 1.0, 0.0, 
  gs = 3.91, -2.678, 
  hw_type = 1
  is_double_j = .true.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 24
  n_block = 0
  n_eigen = 1
  n_restart_vec = 10
  orbs_ratio = 2, 3, 4, 6, 7, 8
&end
EOF
nice ./kshell.exe Cr48_gxpf1a_24.input > log_Cr48_gxpf1a_j24p.txt 2>&1 

rm -f tmp_snapshot_Cr48_gxpf1a_p.ptn_24_* tmp_lv_Cr48_gxpf1a_p.ptn_24_* Cr48_gxpf1a_24.input 


echo "start running log_Cr48_gxpf1a_j28p.txt ..."
cat > Cr48_gxpf1a_28.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5
  fn_int = "gxpf1a.snt"
  fn_ptn = "Cr48_gxpf1a_p.ptn"
  fn_save_wave = "Cr48_gxpf1a_j28p.wav"
  gl = 1.0, 0.0, 
  gs = 3.91, -2.678, 
  hw_type = 1
  is_double_j = .true.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 28
  n_block = 0
  n_eigen = 1
  n_restart_vec = 10
  orbs_ratio = 2, 3, 4, 6, 7, 8
&end
EOF
nice ./kshell.exe Cr48_gxpf1a_28.input > log_Cr48_gxpf1a_j28p.txt 2>&1 

rm -f tmp_snapshot_Cr48_gxpf1a_p.ptn_28_* tmp_lv_Cr48_gxpf1a_p.ptn_28_* Cr48_gxpf1a_28.input 


echo "start running log_Cr48_gxpf1a_j32p.txt ..."
cat > Cr48_gxpf1a_32.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5
  fn_int = "gxpf1a.snt"
  fn_ptn = "Cr48_gxpf1a_p.ptn"
  fn_save_wave = "Cr48_gxpf1a_j32p.wav"
  gl = 1.0, 0.0, 
  gs = 3.91, -2.678, 
  hw_type = 1
  is_double_j = .true.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 32
  n_block = 0
  n_eigen = 1
  n_restart_vec = 10
  orbs_ratio = 2, 3, 4, 6, 7, 8
&end
EOF
nice ./kshell.exe Cr48_gxpf1a_32.input > log_Cr48_gxpf1a_j32p.txt 2>&1 

rm -f tmp_snapshot_Cr48_gxpf1a_p.ptn_32_* tmp_lv_Cr48_gxpf1a_p.ptn_32_* Cr48_gxpf1a_32.input 


./collect_logs.py log_*Cr48_gxpf1a* > summary_Cr48_gxpf1a.txt
echo "Finish computing Cr48_gxpf1a.    See summary_Cr48_gxpf1a.txt"
echo 

