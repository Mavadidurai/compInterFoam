#!/usr/bin/env python3
"""
Diagnostic tool to identify the p_sat calculation bug.
Compares theoretical values vs observed simulation output.
"""

import math

# Parameters from TEST1/system/controlDict phaseChangeCoeffs
Tvap = 2200  # K
L = 9.1e6    # J/kg (latent heat)
R = 174      # J/kg/K (specific gas constant for Ti)
evap_coeff = 0.03
p_ref = 101325  # Pa (assumed standard)

# Observed values from simulation log
T_observed = 8846.3302  # K (max Tl from log)
j_observed = 9.5559451e-09  # kg/m²/s (observed max j_net)

print("="*70)
print(" CLAUSIUS-CLAPEYRON P_SAT DIAGNOSTIC")
print("="*70)
print("\nParameters from controlDict:")
print(f"  Tvap            = {Tvap} K")
print(f"  Latent heat (L) = {L:.2e} J/kg")
print(f"  Gas constant (R)= {R} J/kg/K")
print(f"  Evap coeff (α)  = {evap_coeff}")
print(f"  p_ref (assumed) = {p_ref} Pa")

print(f"\nObserved simulation values:")
print(f"  Temperature T   = {T_observed} K")
print(f"  Mass flux j_net = {j_observed:.4e} kg/m²/s")

# Calculate expected p_sat with CORRECT formula
print("\n" + "="*70)
print(" CORRECT FORMULA: p_sat = p_ref * exp((L/R) * (1/Tvap - 1/T))")
print("="*70)

exponent_correct = (L/R) * (1/Tvap - 1/T_observed)
p_sat_correct = p_ref * math.exp(exponent_correct)

sqrt_term = math.sqrt(2 * math.pi * R * T_observed)
j_expected = evap_coeff * p_sat_correct / sqrt_term

print(f"  Exponent        = {exponent_correct:.4f}")
print(f"  p_sat           = {p_sat_correct:.4e} Pa")
print(f"  sqrt(2πRT)      = {sqrt_term:.2f}")
print(f"  j_evap expected = {j_expected:.4e} kg/m²/s")
print(f"  Ratio (expected/observed) = {j_expected/j_observed:.2e}")

# Calculate what p_sat would be with BACKWARDS formula
print("\n" + "="*70)
print(" BACKWARDS FORMULA: p_sat = p_ref * exp((L/R) * (1/T - 1/Tvap))")
print("="*70)

exponent_backwards = (L/R) * (1/T_observed - 1/Tvap)
p_sat_backwards = p_ref * math.exp(exponent_backwards)
j_backwards = evap_coeff * p_sat_backwards / sqrt_term

print(f"  Exponent        = {exponent_backwards:.4f}")
print(f"  p_sat           = {p_sat_backwards:.4e} Pa")
print(f"  j_evap backwards= {j_backwards:.4e} kg/m²/s")
print(f"  Ratio (backwards/observed) = {j_backwards/j_observed:.2e}")

# Check which matches better
print("\n" + "="*70)
print(" DIAGNOSIS")
print("="*70)

ratio_correct = abs(math.log10(j_expected / j_observed))
ratio_backwards = abs(math.log10(j_backwards / j_observed))

if ratio_backwards < 1:  # Within factor of 10
    print("\n  ⚠️  FOUND THE BUG! ⚠️")
    print("\n  The Clausius-Clapeyron exponent has the WRONG SIGN!")
    print(f"\n  Code is using: (1/T - 1/Tvap) giving exponent = {exponent_backwards:.2f}")
    print(f"  Should be:     (1/Tvap - 1/T) giving exponent = {exponent_correct:.2f}")
    print(f"\n  This causes p_sat to be {p_sat_correct/p_sat_backwards:.2e}× too small!")
    print("\n  FIX: In twoPhaseMixtureThermo.C line ~1140, change:")
    print("       exponent = (L/R)*(inv_Tvap - 1.0/Teval);  // WRONG")
    print("  TO:")
    print("       exponent = (L/R)*(1.0/Teval - inv_Tvap);  // CORRECT")
elif ratio_correct < 1:
    print("\n  The formula appears correct, but simulation doesn't match.")
    print("  Need to investigate other issues (temperature clamping, etc.)")
else:
    print("\n  Neither formula matches observed output.")
    print("  There may be an additional scaling factor or different bug.")

# Work backwards from observed flux to find implied p_sat
print("\n" + "="*70)
print(" REVERSE CALCULATION: What p_sat gives observed j_net?")
print("="*70)

p_sat_implied = (j_observed * sqrt_term) / evap_coeff
exponent_implied = math.log(p_sat_implied / p_ref)
T_inverse_implied = 1/Tvap - exponent_implied * (R/L)

print(f"\n  Implied p_sat            = {p_sat_implied:.4e} Pa")
print(f"  Implied exponent         = {exponent_implied:.4f}")
print(f"  Implied 1/T              = {T_inverse_implied:.8f}")

if T_inverse_implied > 0:
    T_implied = 1/T_inverse_implied
    print(f"  Implied T                = {T_implied:.2f} K")
    if abs(T_implied - T_observed) < 100:
        print("  ⚠️  Temperature matches - likely p_ref or parameter issue")
else:
    print("  ⚠️  Negative 1/T - confirms sign error!")

print("\n" + "="*70)
