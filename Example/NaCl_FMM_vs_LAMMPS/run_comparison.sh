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
FMM_NEAR=$(grep "FMM Direct Coulomb Energy:" fmm_output.log | head -1 | awk '{print $NF}')
FMM_PERIODIC=$(grep "FMM Periodic Correction:" fmm_output.log | head -1 | awk '{print $NF}')

echo "  FMM Direct Coulomb:     $FMM_NEAR kb"
echo "  FMM Periodic Correction: $FMM_PERIODIC kb"
echo "  FMM Total Energy:        $FMM_ENERGY kb"
echo ""

# --- Run LAMMPS ---
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

if [ -n "$FMM_ENERGY" ] && [ -n "$LAMMPS_ENERGY" ]; then
    echo "  FMM    Energy: $FMM_ENERGY kb"
    echo "  LAMMPS Energy: $LAMMPS_ENERGY kb"
    
    # Compute difference using awk
    DIFF=$(echo "$FMM_ENERGY $LAMMPS_ENERGY" | awk '{
        diff = $1 - $2;
        if ($2 != 0) pct = 100.0 * diff / $2;
        else pct = 0;
        printf "  Difference:    %.6f kb (%.4f%%)\n", diff, pct
    }')
    echo "$DIFF"
else
    echo "  Could not extract energies for comparison."
    echo "  Check fmm_output.log and lammps_output.log for errors."
fi

echo "============================================"
echo ""
echo "Full output logs:"
echo "  FMM:    fmm_output.log"
echo "  LAMMPS: lammps_output.log"
