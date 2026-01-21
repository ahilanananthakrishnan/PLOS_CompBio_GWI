READ_ME: GWI INFLAMMATION–HYPERTROPHY–MOTILITY WORKFLOW
=====================================================

This repository implements a three-step finite-element workflow to model
inflammation-driven colon hypertrophy and neuroimmune-mediated motility
changes in Gulf War Illness (GWI). The implementation follows the equations,
notation, and assumptions presented in the accompanying manuscript.

All biological and material parameters are defined in PARAMETERS.inp.

-----------------------------------------------------------------------
OVERVIEW OF MODELING STEPS
-----------------------------------------------------------------------

Step 1: Inflammation-driven hypertrophy
Step 2: Elastic preload (organ-bath baseline stretch)
Step 3: Neural activation under isometric conditions

Each step is executed as a separate Abaqus Standard job and imports the
final configuration from the previous step.

-----------------------------------------------------------------------
STEP 1 — INFLAMMATION-DRIVEN HYPERTROPHY
-----------------------------------------------------------------------
Command:
    abaqus job=STEP1 user=UMAT_STEP1.f

Output:
    STEP1.odb

This hypertrophied configuration is used as the reference state for Step 2.

-----------------------------------------------------------------------
STEP 2 — ELASTIC PRELOAD (ORGAN-BATH BASELINE)
-----------------------------------------------------------------------
Command:
    abaqus job=STEP2 oldjob=STEP1 user=UMAT_STEP2.f

Output:
    STEP2.odb

This prestretched configuration is used as the reference state for Step 3.

-----------------------------------------------------------------------
STEP 3 — NEURAL ACTIVATION (ISOMETRIC)
-----------------------------------------------------------------------
Excitatory activation:
    abaqus job=STEP3 oldjob=STEP2 user=UMAT_STEP3_EXC.f

Inhibitory activation:
    abaqus job=STEP3 oldjob=STEP2 user=UMAT_STEP3_INH.f

-----------------------------------------------------------------------
POST-PROCESSING
-----------------------------------------------------------------------

Active stress output is extracted from Step 3 only.

Command:
    abaqus python active_stress_step3.py

This script computes stress-based metrics from the Step 3 output database
for direct comparison with organ-bath measurements.

-----------------------------------------------------------------------
END
-----------------------------------------------------------------------
