#!/bin/bash
#======================================
# FMM vs LAMMPS Energy Comparison
#
# Runs both FMM and LAMMPS Coulomb calculations
# on the same NaCl configuration and compares
# the initial energies.
#
# Prerequisites:
#   - ClassyMC compiled (standard): for FMM test
#   - ClassyMC compiled with LAMMPS: for LAMMPS test
#     make lammps LAMMPS_LIB_PATH=/path/to/lammps/build
#   - LD_LIBRARY_PATH set if using shared LAMMPS library
#
# Usage:
#   cd Example/NaCl_FMM_vs_LAMMPS
#   bash run_comparison.sh
#======================================

CLASSY="../../classyMC"

# Ensure the KSPACE-enabled LAMMPS library is found
export LD_LIBRARY_PATH=/home/troy/.local/lib:${LD_LIBRARY_PATH:-}

# Check that ClassyMC binary exists
if [ ! -f "$CLASSY" ]; then
    echo "ERROR: ClassyMC binary not found at $CLASSY"
    echo "Please compile ClassyMC first:"
    echo "  make          (for FMM-only test)"
    echo "  make lammps   (for both FMM and LAMMPS tests)"
    exit 1
fi

echo "============================================"
echo " FMM vs LAMMPS Coulomb Energy Comparison"
echo " System: 512 NaCl ions (256 Na+ + 256 Cl-)"
echo " Box: 22.6 x 22.6 x 22.6 Angstroms"
echo " Energy units: kb (Boltzmann)"
echo "============================================"
echo ""

# --- Run FMM ---
echo "--- Running FMM Coulomb calculation ---"
$CLASSY in.nacl_fmm > fmm_output.log 2>&1
FMM_STATUS=$?

if [ $FMM_STATUS -ne 0 ]; then
    echo "WARNING: FMM run exited with status $FMM_STATUS"
fi

# Extract FMM energy
FMM_ENERGY=$(grep "FMM Total Energy:" fmm_output.log | head -1 | awk '{print $NF}')
FMM_NEAR=$(grep "FMM Near-Field Energy:" fmm_output.log | head -1 | awk '{print $NF}')
FMM_FAR=$(grep "FMM Far-Field Energy:" fmm_output.log | head -1 | awk '{print $NF}')
FMM_PERIODIC=$(grep "FMM Periodic Correction:" fmm_output.log | head -1 | awk '{print $NF}')

echo "  FMM Near-Field Energy:   $FMM_NEAR kb"
echo "  FMM Far-Field Energy:    $FMM_FAR kb"
echo "  FMM Periodic Correction: $FMM_PERIODIC kb"
echo "  FMM Total Energy:        $FMM_ENERGY kb"
echo ""

# --- Run Ewald ---
echo "--- Running Ewald Coulomb calculation ---"
$CLASSY in.nacl_ewald > ewald_output.log 2>&1
EWALD_STATUS=$?

if [ $EWALD_STATUS -ne 0 ]; then
    echo "WARNING: Ewald run exited with status $EWALD_STATUS"
fi

# Extract Ewald energy
EWALD_ENERGY=$(grep "Ewald Total Energy:" ewald_output.log | head -1 | awk '{print $NF}')
EWALD_REAL=$(grep "Ewald Real-space Energy:" ewald_output.log | head -1 | awk '{print $NF}')
EWALD_RECIP=$(grep "Ewald Reciprocal Energy:" ewald_output.log | head -1 | awk '{print $NF}')
EWALD_SELF=$(grep "Ewald Self-energy Correction:" ewald_output.log | head -1 | awk '{print $NF}')
EWALD_DIPOLE=$(grep "Ewald Dipole Correction:" ewald_output.log | head -1 | awk '{print $NF}')

echo "  Ewald Real-space:        $EWALD_REAL kb"
echo "  Ewald Reciprocal:        $EWALD_RECIP kb"
echo "  Ewald Self-correction:   $EWALD_SELF kb"
echo "  Ewald Dipole Correction: $EWALD_DIPOLE kb"
echo "  Ewald Total Energy:      $EWALD_ENERGY kb"
echo ""

# --- Run LAMMPS (optional) ---
echo "--- Running LAMMPS Coulomb calculation ---"
$CLASSY in.nacl_lammps > lammps_output.log 2>&1
LAMMPS_STATUS=$?

if [ $LAMMPS_STATUS -ne 0 ]; then
    echo "WARNING: LAMMPS run exited with status $LAMMPS_STATUS"
    echo "  If you see 'Unknown forcefield type', ClassyMC needs to be"
    echo "  compiled with LAMMPS support:"
    echo "    make lammps LAMMPS_LIB_PATH=/path/to/lammps/build"
    echo ""
fi

# Extract LAMMPS energy
LAMMPS_ENERGY=$(grep "Total LAMMPS Energy:" lammps_output.log | head -1 | awk '{print $NF}')

echo "  LAMMPS Total Energy:     $LAMMPS_ENERGY kb"
echo ""

# --- Compare ---
echo "============================================"
echo " COMPARISON"
echo "============================================"

if [ -n "$FMM_ENERGY" ] && [ -n "$EWALD_ENERGY" ]; then
    echo "  FMM   Energy: $FMM_ENERGY kb"
    echo "  Ewald Energy: $EWALD_ENERGY kb"
    
    DIFF=$(echo "$FMM_ENERGY $EWALD_ENERGY" | awk '{
        diff = $1 - $2;
        if ($2 != 0) pct = 100.0 * diff / $2;
        else pct = 0;
        printf "  FMM-Ewald Diff:  %.6f kb (%.4f%%)\n", diff, pct
    }')
    echo "$DIFF"
fi

if [ -n "$FMM_ENERGY" ] && [ -n "$LAMMPS_ENERGY" ]; then
    DIFF=$(echo "$FMM_ENERGY $LAMMPS_ENERGY" | awk '{
        diff = $1 - $2;
        if ($2 != 0) pct = 100.0 * diff / $2;
        else pct = 0;
        printf "  FMM-LAMMPS Diff: %.6f kb (%.4f%%)\n", diff, pct
    }')
    echo "$DIFF"
elif [ -z "$LAMMPS_ENERGY" ]; then
    echo "  (LAMMPS comparison skipped - LAMMPS support not compiled)"
fi

echo "============================================"
echo ""
echo "Full output logs:"
echo "  FMM:    fmm_output.log"
echo "  Ewald:  ewald_output.log"
echo "  LAMMPS: lammps_output.log"
