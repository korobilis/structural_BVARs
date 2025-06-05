# Structural Restrictions in Vector Autoregressions - Replication Code

Replication code for:

**Korobilis, Dimitris (2022), "A new algorithm for structural restrictions in vector autoregressions"**  
*European Economic Review* (forthcoming)  
[SSRN Link](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3557911)

## Structure

### 1. MONTE_CARLO
Replicates Monte Carlo results from Section 3

### 2. EMPIRICAL_PART_A  
Replicates empirical results from Section 4.1 (baseline model)

### 3. EMPIRICAL_PART_B
Replicates empirical results from Section 4.2 (extended model)

## Files

| File | Replicates |
|------|------------|
| `MONTE_CARLO/MONTE_CARLO.m` | Figure 1, Table 2 (Section 3.1), Figures C1-C4 (supplement) |
| `MONTE_CARLO/MONTE_CARLO_TIMES.m` | Table 3 (Section 3.2) |
| `EMPIRICAL_PART_A/BVAR_FSR.m` | Figure 2 (Section 4.1), Figures C6-C9 (supplement) |
| `EMPIRICAL_PART_A/BVAR_FSR_n_r.m` | Figure 3 panel (b) (Section 4.1) |
| `EMPIRICAL_PART_A/BVAR_FSR_1plus4shocks.m` | Figure 3 panel (c) (Section 4.1) |
| `EMPIRICAL_PART_A/BVAR_FSR_1shock.m` | Figure 3 panel (d) (Section 4.1) |
| `EMPIRICAL_PART_A/BVAR_FSR_noninf.m` | Figure C10 (Section C.3.1, supplement) |
| `EMPIRICAL_PART_B/BVAR_FSR.m` | Figure 4 (Section 4.2), Figures C11-C16 (supplement) |

Each MATLAB file contains comments indicating which figures/tables it replicates.
