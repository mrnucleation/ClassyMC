#!/bin/bash
#======================================
# TIP3P FMM vs LAMMPS Energy Comparison
#
# Runs FMM, Ewald, and LAMMPS calculations on the same
# TIP3P water configuration and compares the initial energies.
#
# Prerequisites:
#   - ClassyMC compiled (standard): for FMM and Ewald tests
#   - ClassyMC compiled with LAMMPS: for LAMMPS test
#     make lammps LAMMPS_LIB_PATH=/path/to/lammps/build
#   - LD_LIBRARY_PATH set if using shared LAMMPS library
#
# Usage:
#   cd Example/TIP3P_FMM_vs_LAMMPS
#   bash run_comparison.sh
#======================================

CLASSY="../../classyMC"

# Ensure the KSPACE-enabled LAMMPS library is found
export LD_LIBRARY_PATH=/home/troy/.local/lib:${LD_LIBRARY_PATH:-}

# Check that ClassyMC binary exists
if [ ! -f "$CLASSY" ]; then
    echo "ERROR: ClassyMC binary not found at $CLASSY"
    echo "Please compile ClassyMC first:"
    echo "  make          (for FMM and Ewald tests)"
    echo "  make lammps   (for all three tests)"
    exit 1
fi

echo "============================================"
echo " FMM vs Ewald vs LAMMPS Energy Comparison"
echo " System: 500 TIP3P water molecules"
echo " Box: 25.0 x 25.0 x 25.0 Angstroms"
echo " Energy units: kb (Boltzmann)"
echo "============================================"
echo ""

# --- Run FMM ---
echo "--- Running FMM Hybrid (LJ + FMM) calculation ---"
$CLASSY in.tip3p_fmm > fmm_output.log 2>&1
FMM_STATUS=$?

if [ $FMM_STATUS -ne 0 ]; then
    echo "WARNING: FMM run exited with status $FMM_STATUS"
fi

# Extract FMM energy components
FMM_LJ=$(grep "Total Pair Energy:" fmm_output.log | head -1 | awk '{print $NF}')
FMM_TAIL=$(grep "Total Tail Corrections:" fmm_output.log | head -1 | awk '{print $NF}')
FMM_NEAR=$(grep "FMM Near-Field Energy:" fmm_output.log | head -1 | awk '{print $NF}')
FMM_FAR=$(grep "FMM Far-Field Energy:" fmm_output.log | head -1 | awk '{print $NF}')
FMM_LATTICE=$(grep "FMM Periodic Lattice-Image Correction" fmm_output.log | head -1 | awk '{print $NF}')
FMM_SURFACE=$(grep "FMM Surface Dipole Correction" fmm_output.log | head -1 | awk '{print $NF}')
FMM_COULOMB=$(grep "FMM Total Energy:" fmm_output.log | head -1 | awk '{print $NF}')
FMM_TOTAL=$(grep "Hybrid Forcefield Energy:" fmm_output.log | head -1 | awk '{print $NF}')

echo "  LJ Pair Energy:          $FMM_LJ kb"
echo "  LJ Tail Corrections:     $FMM_TAIL kb"
echo "  FMM Near-Field Energy:   $FMM_NEAR kb"
echo "  FMM Far-Field Energy:    $FMM_FAR kb"
echo "  FMM Lattice-Image Corr:  $FMM_LATTICE kb"
echo "  FMM Surface Dipole:      $FMM_SURFACE kb"
echo "  FMM Coulomb Total:       $FMM_COULOMB kb"
echo "  FMM Hybrid Total:        $FMM_TOTAL kb"
echo ""

# --- Run Ewald ---
echo "--- Running LJ + Ewald calculation ---"
$CLASSY in.tip3p_ewald > ewald_output.log 2>&1
EWALD_STATUS=$?

if [ $EWALD_STATUS -ne 0 ]; then
    echo "WARNING: Ewald run exited with status $EWALD_STATUS"
fi

# Extract Ewald energy components
EWALD_LJ=$(grep "LJ Energy:" ewald_output.log | head -1 | awk '{print $NF}')
EWALD_REAL=$(grep "Coulomb Real-space:" ewald_output.log | head -1 | awk '{print $NF}')
EWALD_RECIP=$(grep "Coulomb Reciprocal:" ewald_output.log | head -1 | awk '{print $NF}')
EWALD_SELF=$(grep "Coulomb Self-correction:" ewald_output.log | head -1 | awk '{print $NF}')
EWALD_DIPOLE=$(grep "Coulomb Surface Dipole Correction" ewald_output.log | head -1 | awk '{print $NF}')
EWALD_EXCL=$(grep "Coulomb Exclusion Correction:" ewald_output.log | head -1 | awk '{print $NF}')
EWALD_TOTAL=$(grep "Total Energy:" ewald_output.log | head -1 | awk '{print $NF}')

echo "  LJ Energy:               $EWALD_LJ kb"
echo "  Coulomb Real-space:      $EWALD_REAL kb"
echo "  Coulomb Reciprocal:      $EWALD_RECIP kb"
echo "  Coulomb Self-correction: $EWALD_SELF kb"
echo "  Coulomb Dipole Corr:     $EWALD_DIPOLE kb"
echo "  Coulomb Exclusion Corr:  $EWALD_EXCL kb"
echo "  Ewald Total Energy:      $EWALD_TOTAL kb"
echo ""

# --- Run LAMMPS ---
echo "--- Running LAMMPS (lj/cut/coul/long + PPPM) calculation ---"
$CLASSY in.tip3p_lammps > lammps_output.log 2>&1
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
echo " COMPARISON (Hybrid Total Energies)"
echo "============================================"

if [ -n "$FMM_TOTAL" ] && [ -n "$EWALD_TOTAL" ]; then
    echo "  FMM   Total: $FMM_TOTAL kb"
    echo "  Ewald Total: $EWALD_TOTAL kb"

    DIFF=$(echo "$FMM_TOTAL $EWALD_TOTAL" | awk '{
        diff = $1 - $2;
        if ($2 != 0) pct = 100.0 * diff / $2;
        else pct = 0;
        printf "  FMM-Ewald Diff:  %.6f kb (%.4f%%)\n", diff, pct
    }')
    echo "$DIFF"
fi

if [ -n "$FMM_TOTAL" ] && [ -n "$LAMMPS_ENERGY" ]; then
    DIFF=$(echo "$FMM_TOTAL $LAMMPS_ENERGY" | awk '{
        diff = $1 - $2;
        if ($2 != 0) pct = 100.0 * diff / $2;
        else pct = 0;
        printf "  FMM-LAMMPS Diff: %.6f kb (%.4f%%)\n", diff, pct
    }')
    echo "$DIFF"
elif [ -z "$LAMMPS_ENERGY" ]; then
    echo "  (LAMMPS comparison skipped - LAMMPS support not compiled)"
fi

if [ -n "$EWALD_TOTAL" ] && [ -n "$LAMMPS_ENERGY" ]; then
    DIFF=$(echo "$EWALD_TOTAL $LAMMPS_ENERGY" | awk '{
        diff = $1 - $2;
        if ($2 != 0) pct = 100.0 * diff / $2;
        else pct = 0;
        printf "  Ewald-LAMMPS Diff: %.6f kb (%.4f%%)\n", diff, pct
    }')
    echo "$DIFF"
fi

echo ""
echo "Notes:"
echo "  - All three methods use tin-foil (eps_s=inf) boundary conditions (ClassyMC"
echo "    default 'surface tinfoil', and LAMMPS PPPM default). No surface dipole term."
echo "  - LAMMPS is run with atom_style charge (no bond info available), so its raw"
echo "    output contains the full Coulomb of all 1500 charges -- including intramolecular"
echo "    (see 'Raw LAMMPS Energy' and 'LAMMPS Intra-Coulomb (subtracted)' in lammps_output.log)."
echo "    ClassyMC subtracts the intra term so LAMMPS reports the same inter-molecular only"
echo "    energy that Ewald and FMM compute."
echo "  - FMM uses 'image_shells N' to sum the Coulomb contribution from N periodic image"
echo "    shells. Cubic shell summation converges to the vacuum (eps_s=1) BC; ClassyMC"
echo "    strips the surface dipole 2*pi*|M|^2/(3V) automatically to report the tin-foil total."
echo "  - Residual FMM-Ewald discrepancy ~0.003% is truncation error of 3 image shells."
echo "  - Residual Ewald-LAMMPS ~0.17% is PPPM grid/order precision."

echo "============================================"
echo ""
echo "Full output logs:"
echo "  FMM:    fmm_output.log"
echo "  Ewald:  ewald_output.log"
echo "  LAMMPS: lammps_output.log"
