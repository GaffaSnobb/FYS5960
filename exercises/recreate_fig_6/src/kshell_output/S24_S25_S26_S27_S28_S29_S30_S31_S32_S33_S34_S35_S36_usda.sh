#!/bin/sh 
# export OMP_STACKSIZE=1g
export GFORTRAN_UNBUFFERED_PRECONNECTED=y
# ulimit -s unlimited

# ---------- S24_usda --------------
echo "start running log_S24_usda_m0p.txt ..."
cat > S24_usda_0.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5, 
  fn_int = "usda.snt"
  fn_ptn = "S24_usda_p.ptn"
  fn_save_wave = "S24_usda_m0p.wav"
  gl = 1.0, 0.0, 
  gs = 3.91, -2.678, 
  hw_type = 2
  is_double_j = .false.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 0
  n_block = 0
  n_eigen = 100
  n_restart_vec = 150
&end
EOF
nice ./kshell.exe S24_usda_0.input > log_S24_usda_m0p.txt 2>&1 

rm -f tmp_snapshot_S24_usda_p.ptn_0_* tmp_lv_S24_usda_p.ptn_0_* S24_usda_0.input 


# --------------- transition probabilities --------------

echo "start running log_S24_usda_tr_m0p_m0p.txt ..."
cat > S24_usda_0_0.input <<EOF
&input
  fn_int   = "usda.snt"
  fn_ptn_l = "S24_usda_p.ptn"
  fn_ptn_r = "S24_usda_p.ptn"
  fn_load_wave_l = "S24_usda_m0p.wav"
  fn_load_wave_r = "S24_usda_m0p.wav"
  hw_type = 2
  eff_charge = 1.5, 0.5
  gl = 1.0, 0.0
  gs = 3.91, -2.678
&end
EOF
nice ./transit.exe S24_usda_0_0.input > log_S24_usda_tr_m0p_m0p.txt 2>&1 

rm -f S24_usda_0_0.input


./collect_logs.py log_*S24_usda* > summary_S24_usda.txt
echo "Finish computing S24_usda.    See summary_S24_usda.txt"
echo 

# ---------- S25_usda --------------
echo "start running log_S25_usda_m1p.txt ..."
cat > S25_usda_1.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5, 
  fn_int = "usda.snt"
  fn_ptn = "S25_usda_p.ptn"
  fn_save_wave = "S25_usda_m1p.wav"
  gl = 1.0, 0.0, 
  gs = 3.91, -2.678, 
  hw_type = 2
  is_double_j = .false.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 1
  n_block = 0
  n_eigen = 100
  n_restart_vec = 150
&end
EOF
nice ./kshell.exe S25_usda_1.input > log_S25_usda_m1p.txt 2>&1 

rm -f tmp_snapshot_S25_usda_p.ptn_1_* tmp_lv_S25_usda_p.ptn_1_* S25_usda_1.input 


# --------------- transition probabilities --------------

echo "start running log_S25_usda_tr_m1p_m1p.txt ..."
cat > S25_usda_1_1.input <<EOF
&input
  fn_int   = "usda.snt"
  fn_ptn_l = "S25_usda_p.ptn"
  fn_ptn_r = "S25_usda_p.ptn"
  fn_load_wave_l = "S25_usda_m1p.wav"
  fn_load_wave_r = "S25_usda_m1p.wav"
  hw_type = 2
  eff_charge = 1.5, 0.5
  gl = 1.0, 0.0
  gs = 3.91, -2.678
&end
EOF
nice ./transit.exe S25_usda_1_1.input > log_S25_usda_tr_m1p_m1p.txt 2>&1 

rm -f S25_usda_1_1.input


./collect_logs.py log_*S25_usda* > summary_S25_usda.txt
echo "Finish computing S25_usda.    See summary_S25_usda.txt"
echo 

# ---------- S26_usda --------------
echo "start running log_S26_usda_m0p.txt ..."
cat > S26_usda_0.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5, 
  fn_int = "usda.snt"
  fn_ptn = "S26_usda_p.ptn"
  fn_save_wave = "S26_usda_m0p.wav"
  gl = 1.0, 0.0, 
  gs = 3.91, -2.678, 
  hw_type = 2
  is_double_j = .false.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 0
  n_block = 0
  n_eigen = 100
  n_restart_vec = 150
&end
EOF
nice ./kshell.exe S26_usda_0.input > log_S26_usda_m0p.txt 2>&1 

rm -f tmp_snapshot_S26_usda_p.ptn_0_* tmp_lv_S26_usda_p.ptn_0_* S26_usda_0.input 


# --------------- transition probabilities --------------

echo "start running log_S26_usda_tr_m0p_m0p.txt ..."
cat > S26_usda_0_0.input <<EOF
&input
  fn_int   = "usda.snt"
  fn_ptn_l = "S26_usda_p.ptn"
  fn_ptn_r = "S26_usda_p.ptn"
  fn_load_wave_l = "S26_usda_m0p.wav"
  fn_load_wave_r = "S26_usda_m0p.wav"
  hw_type = 2
  eff_charge = 1.5, 0.5
  gl = 1.0, 0.0
  gs = 3.91, -2.678
&end
EOF
nice ./transit.exe S26_usda_0_0.input > log_S26_usda_tr_m0p_m0p.txt 2>&1 

rm -f S26_usda_0_0.input


./collect_logs.py log_*S26_usda* > summary_S26_usda.txt
echo "Finish computing S26_usda.    See summary_S26_usda.txt"
echo 

# ---------- S27_usda --------------
echo "start running log_S27_usda_m1p.txt ..."
cat > S27_usda_1.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5, 
  fn_int = "usda.snt"
  fn_ptn = "S27_usda_p.ptn"
  fn_save_wave = "S27_usda_m1p.wav"
  gl = 1.0, 0.0, 
  gs = 3.91, -2.678, 
  hw_type = 2
  is_double_j = .false.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 1
  n_block = 0
  n_eigen = 100
  n_restart_vec = 150
&end
EOF
nice ./kshell.exe S27_usda_1.input > log_S27_usda_m1p.txt 2>&1 

rm -f tmp_snapshot_S27_usda_p.ptn_1_* tmp_lv_S27_usda_p.ptn_1_* S27_usda_1.input 


# --------------- transition probabilities --------------

echo "start running log_S27_usda_tr_m1p_m1p.txt ..."
cat > S27_usda_1_1.input <<EOF
&input
  fn_int   = "usda.snt"
  fn_ptn_l = "S27_usda_p.ptn"
  fn_ptn_r = "S27_usda_p.ptn"
  fn_load_wave_l = "S27_usda_m1p.wav"
  fn_load_wave_r = "S27_usda_m1p.wav"
  hw_type = 2
  eff_charge = 1.5, 0.5
  gl = 1.0, 0.0
  gs = 3.91, -2.678
&end
EOF
nice ./transit.exe S27_usda_1_1.input > log_S27_usda_tr_m1p_m1p.txt 2>&1 

rm -f S27_usda_1_1.input


./collect_logs.py log_*S27_usda* > summary_S27_usda.txt
echo "Finish computing S27_usda.    See summary_S27_usda.txt"
echo 

# ---------- S28_usda --------------
echo "start running log_S28_usda_m0p.txt ..."
cat > S28_usda_0.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5, 
  fn_int = "usda.snt"
  fn_ptn = "S28_usda_p.ptn"
  fn_save_wave = "S28_usda_m0p.wav"
  gl = 1.0, 0.0, 
  gs = 3.91, -2.678, 
  hw_type = 2
  is_double_j = .false.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 0
  n_block = 0
  n_eigen = 100
  n_restart_vec = 150
&end
EOF
nice ./kshell.exe S28_usda_0.input > log_S28_usda_m0p.txt 2>&1 

rm -f tmp_snapshot_S28_usda_p.ptn_0_* tmp_lv_S28_usda_p.ptn_0_* S28_usda_0.input 


# --------------- transition probabilities --------------

echo "start running log_S28_usda_tr_m0p_m0p.txt ..."
cat > S28_usda_0_0.input <<EOF
&input
  fn_int   = "usda.snt"
  fn_ptn_l = "S28_usda_p.ptn"
  fn_ptn_r = "S28_usda_p.ptn"
  fn_load_wave_l = "S28_usda_m0p.wav"
  fn_load_wave_r = "S28_usda_m0p.wav"
  hw_type = 2
  eff_charge = 1.5, 0.5
  gl = 1.0, 0.0
  gs = 3.91, -2.678
&end
EOF
nice ./transit.exe S28_usda_0_0.input > log_S28_usda_tr_m0p_m0p.txt 2>&1 

rm -f S28_usda_0_0.input


./collect_logs.py log_*S28_usda* > summary_S28_usda.txt
echo "Finish computing S28_usda.    See summary_S28_usda.txt"
echo 

# ---------- S29_usda --------------
echo "start running log_S29_usda_m1p.txt ..."
cat > S29_usda_1.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5, 
  fn_int = "usda.snt"
  fn_ptn = "S29_usda_p.ptn"
  fn_save_wave = "S29_usda_m1p.wav"
  gl = 1.0, 0.0, 
  gs = 3.91, -2.678, 
  hw_type = 2
  is_double_j = .false.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 1
  n_block = 0
  n_eigen = 100
  n_restart_vec = 150
&end
EOF
nice ./kshell.exe S29_usda_1.input > log_S29_usda_m1p.txt 2>&1 

rm -f tmp_snapshot_S29_usda_p.ptn_1_* tmp_lv_S29_usda_p.ptn_1_* S29_usda_1.input 


# --------------- transition probabilities --------------

echo "start running log_S29_usda_tr_m1p_m1p.txt ..."
cat > S29_usda_1_1.input <<EOF
&input
  fn_int   = "usda.snt"
  fn_ptn_l = "S29_usda_p.ptn"
  fn_ptn_r = "S29_usda_p.ptn"
  fn_load_wave_l = "S29_usda_m1p.wav"
  fn_load_wave_r = "S29_usda_m1p.wav"
  hw_type = 2
  eff_charge = 1.5, 0.5
  gl = 1.0, 0.0
  gs = 3.91, -2.678
&end
EOF
nice ./transit.exe S29_usda_1_1.input > log_S29_usda_tr_m1p_m1p.txt 2>&1 

rm -f S29_usda_1_1.input


./collect_logs.py log_*S29_usda* > summary_S29_usda.txt
echo "Finish computing S29_usda.    See summary_S29_usda.txt"
echo 

# ---------- S30_usda --------------
echo "start running log_S30_usda_m0p.txt ..."
cat > S30_usda_0.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5, 
  fn_int = "usda.snt"
  fn_ptn = "S30_usda_p.ptn"
  fn_save_wave = "S30_usda_m0p.wav"
  gl = 1.0, 0.0, 
  gs = 3.91, -2.678, 
  hw_type = 2
  is_double_j = .false.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 0
  n_block = 0
  n_eigen = 100
  n_restart_vec = 150
&end
EOF
nice ./kshell.exe S30_usda_0.input > log_S30_usda_m0p.txt 2>&1 

rm -f tmp_snapshot_S30_usda_p.ptn_0_* tmp_lv_S30_usda_p.ptn_0_* S30_usda_0.input 


# --------------- transition probabilities --------------

echo "start running log_S30_usda_tr_m0p_m0p.txt ..."
cat > S30_usda_0_0.input <<EOF
&input
  fn_int   = "usda.snt"
  fn_ptn_l = "S30_usda_p.ptn"
  fn_ptn_r = "S30_usda_p.ptn"
  fn_load_wave_l = "S30_usda_m0p.wav"
  fn_load_wave_r = "S30_usda_m0p.wav"
  hw_type = 2
  eff_charge = 1.5, 0.5
  gl = 1.0, 0.0
  gs = 3.91, -2.678
&end
EOF
nice ./transit.exe S30_usda_0_0.input > log_S30_usda_tr_m0p_m0p.txt 2>&1 

rm -f S30_usda_0_0.input


./collect_logs.py log_*S30_usda* > summary_S30_usda.txt
echo "Finish computing S30_usda.    See summary_S30_usda.txt"
echo 

# ---------- S31_usda --------------
echo "start running log_S31_usda_m1p.txt ..."
cat > S31_usda_1.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5, 
  fn_int = "usda.snt"
  fn_ptn = "S31_usda_p.ptn"
  fn_save_wave = "S31_usda_m1p.wav"
  gl = 1.0, 0.0, 
  gs = 3.91, -2.678, 
  hw_type = 2
  is_double_j = .false.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 1
  n_block = 0
  n_eigen = 100
  n_restart_vec = 150
&end
EOF
nice ./kshell.exe S31_usda_1.input > log_S31_usda_m1p.txt 2>&1 

rm -f tmp_snapshot_S31_usda_p.ptn_1_* tmp_lv_S31_usda_p.ptn_1_* S31_usda_1.input 


# --------------- transition probabilities --------------

echo "start running log_S31_usda_tr_m1p_m1p.txt ..."
cat > S31_usda_1_1.input <<EOF
&input
  fn_int   = "usda.snt"
  fn_ptn_l = "S31_usda_p.ptn"
  fn_ptn_r = "S31_usda_p.ptn"
  fn_load_wave_l = "S31_usda_m1p.wav"
  fn_load_wave_r = "S31_usda_m1p.wav"
  hw_type = 2
  eff_charge = 1.5, 0.5
  gl = 1.0, 0.0
  gs = 3.91, -2.678
&end
EOF
nice ./transit.exe S31_usda_1_1.input > log_S31_usda_tr_m1p_m1p.txt 2>&1 

rm -f S31_usda_1_1.input


./collect_logs.py log_*S31_usda* > summary_S31_usda.txt
echo "Finish computing S31_usda.    See summary_S31_usda.txt"
echo 

# ---------- S32_usda --------------
echo "start running log_S32_usda_m0p.txt ..."
cat > S32_usda_0.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5, 
  fn_int = "usda.snt"
  fn_ptn = "S32_usda_p.ptn"
  fn_save_wave = "S32_usda_m0p.wav"
  gl = 1.0, 0.0, 
  gs = 3.91, -2.678, 
  hw_type = 2
  is_double_j = .false.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 0
  n_block = 0
  n_eigen = 100
  n_restart_vec = 150
&end
EOF
nice ./kshell.exe S32_usda_0.input > log_S32_usda_m0p.txt 2>&1 

rm -f tmp_snapshot_S32_usda_p.ptn_0_* tmp_lv_S32_usda_p.ptn_0_* S32_usda_0.input 


# --------------- transition probabilities --------------

echo "start running log_S32_usda_tr_m0p_m0p.txt ..."
cat > S32_usda_0_0.input <<EOF
&input
  fn_int   = "usda.snt"
  fn_ptn_l = "S32_usda_p.ptn"
  fn_ptn_r = "S32_usda_p.ptn"
  fn_load_wave_l = "S32_usda_m0p.wav"
  fn_load_wave_r = "S32_usda_m0p.wav"
  hw_type = 2
  eff_charge = 1.5, 0.5
  gl = 1.0, 0.0
  gs = 3.91, -2.678
&end
EOF
nice ./transit.exe S32_usda_0_0.input > log_S32_usda_tr_m0p_m0p.txt 2>&1 

rm -f S32_usda_0_0.input


./collect_logs.py log_*S32_usda* > summary_S32_usda.txt
echo "Finish computing S32_usda.    See summary_S32_usda.txt"
echo 

# ---------- S33_usda --------------
echo "start running log_S33_usda_m1p.txt ..."
cat > S33_usda_1.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5, 
  fn_int = "usda.snt"
  fn_ptn = "S33_usda_p.ptn"
  fn_save_wave = "S33_usda_m1p.wav"
  gl = 1.0, 0.0, 
  gs = 3.91, -2.678, 
  hw_type = 2
  is_double_j = .false.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 1
  n_block = 0
  n_eigen = 100
  n_restart_vec = 150
&end
EOF
nice ./kshell.exe S33_usda_1.input > log_S33_usda_m1p.txt 2>&1 

rm -f tmp_snapshot_S33_usda_p.ptn_1_* tmp_lv_S33_usda_p.ptn_1_* S33_usda_1.input 


# --------------- transition probabilities --------------

echo "start running log_S33_usda_tr_m1p_m1p.txt ..."
cat > S33_usda_1_1.input <<EOF
&input
  fn_int   = "usda.snt"
  fn_ptn_l = "S33_usda_p.ptn"
  fn_ptn_r = "S33_usda_p.ptn"
  fn_load_wave_l = "S33_usda_m1p.wav"
  fn_load_wave_r = "S33_usda_m1p.wav"
  hw_type = 2
  eff_charge = 1.5, 0.5
  gl = 1.0, 0.0
  gs = 3.91, -2.678
&end
EOF
nice ./transit.exe S33_usda_1_1.input > log_S33_usda_tr_m1p_m1p.txt 2>&1 

rm -f S33_usda_1_1.input


./collect_logs.py log_*S33_usda* > summary_S33_usda.txt
echo "Finish computing S33_usda.    See summary_S33_usda.txt"
echo 

# ---------- S34_usda --------------
echo "start running log_S34_usda_m0p.txt ..."
cat > S34_usda_0.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5, 
  fn_int = "usda.snt"
  fn_ptn = "S34_usda_p.ptn"
  fn_save_wave = "S34_usda_m0p.wav"
  gl = 1.0, 0.0, 
  gs = 3.91, -2.678, 
  hw_type = 2
  is_double_j = .false.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 0
  n_block = 0
  n_eigen = 100
  n_restart_vec = 150
&end
EOF
nice ./kshell.exe S34_usda_0.input > log_S34_usda_m0p.txt 2>&1 

rm -f tmp_snapshot_S34_usda_p.ptn_0_* tmp_lv_S34_usda_p.ptn_0_* S34_usda_0.input 


# --------------- transition probabilities --------------

echo "start running log_S34_usda_tr_m0p_m0p.txt ..."
cat > S34_usda_0_0.input <<EOF
&input
  fn_int   = "usda.snt"
  fn_ptn_l = "S34_usda_p.ptn"
  fn_ptn_r = "S34_usda_p.ptn"
  fn_load_wave_l = "S34_usda_m0p.wav"
  fn_load_wave_r = "S34_usda_m0p.wav"
  hw_type = 2
  eff_charge = 1.5, 0.5
  gl = 1.0, 0.0
  gs = 3.91, -2.678
&end
EOF
nice ./transit.exe S34_usda_0_0.input > log_S34_usda_tr_m0p_m0p.txt 2>&1 

rm -f S34_usda_0_0.input


./collect_logs.py log_*S34_usda* > summary_S34_usda.txt
echo "Finish computing S34_usda.    See summary_S34_usda.txt"
echo 

# ---------- S35_usda --------------
echo "start running log_S35_usda_m1p.txt ..."
cat > S35_usda_1.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5, 
  fn_int = "usda.snt"
  fn_ptn = "S35_usda_p.ptn"
  fn_save_wave = "S35_usda_m1p.wav"
  gl = 1.0, 0.0, 
  gs = 3.91, -2.678, 
  hw_type = 2
  is_double_j = .false.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 1
  n_block = 0
  n_eigen = 100
  n_restart_vec = 150
&end
EOF
nice ./kshell.exe S35_usda_1.input > log_S35_usda_m1p.txt 2>&1 

rm -f tmp_snapshot_S35_usda_p.ptn_1_* tmp_lv_S35_usda_p.ptn_1_* S35_usda_1.input 


# --------------- transition probabilities --------------

echo "start running log_S35_usda_tr_m1p_m1p.txt ..."
cat > S35_usda_1_1.input <<EOF
&input
  fn_int   = "usda.snt"
  fn_ptn_l = "S35_usda_p.ptn"
  fn_ptn_r = "S35_usda_p.ptn"
  fn_load_wave_l = "S35_usda_m1p.wav"
  fn_load_wave_r = "S35_usda_m1p.wav"
  hw_type = 2
  eff_charge = 1.5, 0.5
  gl = 1.0, 0.0
  gs = 3.91, -2.678
&end
EOF
nice ./transit.exe S35_usda_1_1.input > log_S35_usda_tr_m1p_m1p.txt 2>&1 

rm -f S35_usda_1_1.input


./collect_logs.py log_*S35_usda* > summary_S35_usda.txt
echo "Finish computing S35_usda.    See summary_S35_usda.txt"
echo 

# ---------- S36_usda --------------
echo "start running log_S36_usda_m0p.txt ..."
cat > S36_usda_0.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5, 
  fn_int = "usda.snt"
  fn_ptn = "S36_usda_p.ptn"
  fn_save_wave = "S36_usda_m0p.wav"
  gl = 1.0, 0.0, 
  gs = 3.91, -2.678, 
  hw_type = 2
  is_double_j = .false.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 0
  n_block = 0
  n_eigen = 100
  n_restart_vec = 150
&end
EOF
nice ./kshell.exe S36_usda_0.input > log_S36_usda_m0p.txt 2>&1 

rm -f tmp_snapshot_S36_usda_p.ptn_0_* tmp_lv_S36_usda_p.ptn_0_* S36_usda_0.input 


# --------------- transition probabilities --------------

echo "start running log_S36_usda_tr_m0p_m0p.txt ..."
cat > S36_usda_0_0.input <<EOF
&input
  fn_int   = "usda.snt"
  fn_ptn_l = "S36_usda_p.ptn"
  fn_ptn_r = "S36_usda_p.ptn"
  fn_load_wave_l = "S36_usda_m0p.wav"
  fn_load_wave_r = "S36_usda_m0p.wav"
  hw_type = 2
  eff_charge = 1.5, 0.5
  gl = 1.0, 0.0
  gs = 3.91, -2.678
&end
EOF
nice ./transit.exe S36_usda_0_0.input > log_S36_usda_tr_m0p_m0p.txt 2>&1 

rm -f S36_usda_0_0.input


./collect_logs.py log_*S36_usda* > summary_S36_usda.txt
echo "Finish computing S36_usda.    See summary_S36_usda.txt"
echo 

