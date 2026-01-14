#!/bin/bash
#SBATCH --job-name=concatenate_full_v2
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=1:30:00
#SBATCH --error=../slurm_log/%x_%A.err
#SBATCH --output=../slurm_log/%x_%A.out

# Load environment
module purge
module load anaconda3/2024.10
conda activate scseq


# bash ./concatenate.sh PA14+BC:A4,A5,A6,B4,B5,B6,C4,C5,C6,D4,D5,D6,E4,E5,E6,F4,F5,F6,G4,G5,G6,H4,H5,H6 ../rawdata/r1r3_welldemultiplex/r3r1/ 31 ../rawdata/condition_demultiplex/ &
# bash ./concatenate.sh PA14+SA:A1,A2,A3,B1,B2,B3,C1,C2,C3,D1,D2,D3,E1,E2,E3,F1,F2,F3,G1,G2,G3,H1,H2,H3 ../rawdata/r1r3_welldemultiplex/r3r1/ 31 ../rawdata/condition_demultiplex/ &
# bash ./concatenate.sh PAO1+BC:A10,A11,A12,B10,B11,B12,C10,C11,C12,D10,D11,D12,E10,E11,E12,F10,F11,F12,G10,G11,G12,H10,H11,H12 ../rawdata/r1r3_welldemultiplex/r3r1/ 31 ../rawdata/condition_demultiplex/ &
# bash ./concatenate.sh PAO1+SA:A7,A8,A9,B7,B8,B9,C7,C8,C9,D7,D8,D9,E7,E8,E9,F7,F8,F9,G7,G8,G9,H7,H8,H9 ../rawdata/r1r3_welldemultiplex/r3r1/ 31 ../rawdata/condition_demultiplex/ &

# bash ./concatenate.sh PA14+BC:A4,A5,A6,B4,B5,B6,C4,C5,C6,D4,D5,D6,E4,E5,E6,F4,F5,F6,G4,G5,G6,H4,H5,H6 ../rawdata/r1r4_welldemultiplex/r4/ 4 ../rawdata/condition_demultiplex/ &
# bash ./concatenate.sh PA14+SA:A1,A2,A3,B1,B2,B3,C1,C2,C3,D1,D2,D3,E1,E2,E3,F1,F2,F3,G1,G2,G3,H1,H2,H3 ../rawdata/r1r4_welldemultiplex/r4/ 4 ../rawdata/condition_demultiplex/ &
# bash ./concatenate.sh PAO1+BC:A10,A11,A12,B10,B11,B12,C10,C11,C12,D10,D11,D12,E10,E11,E12,F10,F11,F12,G10,G11,G12,H10,H11,H12 ../rawdata/r1r4_welldemultiplex/r4/ 4 ../rawdata/condition_demultiplex/ &
# bash ./concatenate.sh PAO1+SA:A7,A8,A9,B7,B8,B9,C7,C8,C9,D7,D8,D9,E7,E8,E9,F7,F8,F9,G7,G8,G9,H7,H8,H9 ../rawdata/r1r4_welldemultiplex/r4/ 4 ../rawdata/condition_demultiplex/ &

# bash ./concatenate.sh PA14+BC:A4,A5,A6,B4,B5,B6,C4,C5,C6,D4,D5,D6,E4,E5,E6,F4,F5,F6,G4,G5,G6,H4,H5,H6 ../rawdata/r1r3_welldemultiplex_full/r3r1/ 31 ../rawdata/condition_demultiplex_full/ &
# bash ./concatenate.sh PA14+SA:A1,A2,A3,B1,B2,B3,C1,C2,C3,D1,D2,D3,E1,E2,E3,F1,F2,F3,G1,G2,G3,H1,H2,H3 ../rawdata/r1r3_welldemultiplex_full/r3r1/ 31 ../rawdata/condition_demultiplex_full/ &
# bash ./concatenate.sh PAO1+BC:A10,A11,A12,B10,B11,B12,C10,C11,C12,D10,D11,D12,E10,E11,E12,F10,F11,F12,G10,G11,G12,H10,H11,H12 ../rawdata/r1r3_welldemultiplex_full/r3r1/ 31 ../rawdata/condition_demultiplex_full/ &
# bash ./concatenate.sh PAO1+SA:A7,A8,A9,B7,B8,B9,C7,C8,C9,D7,D8,D9,E7,E8,E9,F7,F8,F9,G7,G8,G9,H7,H8,H9 ../rawdata/r1r3_welldemultiplex_full/r3r1/ 31 ../rawdata/condition_demultiplex_full/ &

# bash ./concatenate.sh PA14+BC:A4,A5,A6,B4,B5,B6,C4,C5,C6,D4,D5,D6,E4,E5,E6,F4,F5,F6,G4,G5,G6,H4,H5,H6 ../rawdata/r1r4_welldemultiplex_full/r4/ 4 ../rawdata/condition_demultiplex_full/ &
# bash ./concatenate.sh PA14+SA:A1,A2,A3,B1,B2,B3,C1,C2,C3,D1,D2,D3,E1,E2,E3,F1,F2,F3,G1,G2,G3,H1,H2,H3 ../rawdata/r1r4_welldemultiplex_full/r4/ 4 ../rawdata/condition_demultiplex_full/ &
# bash ./concatenate.sh PAO1+BC:A10,A11,A12,B10,B11,B12,C10,C11,C12,D10,D11,D12,E10,E11,E12,F10,F11,F12,G10,G11,G12,H10,H11,H12 ../rawdata/r1r4_welldemultiplex_full/r4/ 4 ../rawdata/condition_demultiplex_full/ &
# bash ./concatenate.sh PAO1+SA:A7,A8,A9,B7,B8,B9,C7,C8,C9,D7,D8,D9,E7,E8,E9,F7,F8,F9,G7,G8,G9,H7,H8,H9 ../rawdata/r1r4_welldemultiplex_full/r4/ 4 ../rawdata/condition_demultiplex_full/ &

bash ./concatenate.sh PA14+BC:A4,A5,A6,B4,B5,B6,C4,C5,C6,D4,D5,D6,E4,E5,E6,F4,F5,F6,G4,G5,G6,H4,H5,H6 ../rawdata/r1r3_welldemultiplex_full_v2/r3r1/ 31 ../rawdata/condition_demultiplex_full_v2/ &
bash ./concatenate.sh PA14+SA:A1,A2,A3,B1,B2,B3,C1,C2,C3,D1,D2,D3,E1,E2,E3,F1,F2,F3,G1,G2,G3,H1,H2,H3 ../rawdata/r1r3_welldemultiplex_full_v2/r3r1/ 31 ../rawdata/condition_demultiplex_full_v2/ &
bash ./concatenate.sh PAO1+BC:A10,A11,A12,B10,B11,B12,C10,C11,C12,D10,D11,D12,E10,E11,E12,F10,F11,F12,G10,G11,G12,H10,H11,H12 ../rawdata/r1r3_welldemultiplex_full_v2/r3r1/ 31 ../rawdata/condition_demultiplex_full_v2/ &
bash ./concatenate.sh PAO1+SA:A7,A8,A9,B7,B8,B9,C7,C8,C9,D7,D8,D9,E7,E8,E9,F7,F8,F9,G7,G8,G9,H7,H8,H9 ../rawdata/r1r3_welldemultiplex_full_v2/r3r1/ 31 ../rawdata/condition_demultiplex_full_v2/ &

bash ./concatenate.sh PA14+BC:A4,A5,A6,B4,B5,B6,C4,C5,C6,D4,D5,D6,E4,E5,E6,F4,F5,F6,G4,G5,G6,H4,H5,H6 ../rawdata/r1r4_welldemultiplex_full_v2/r4/ 4 ../rawdata/condition_demultiplex_full_v2/ &
bash ./concatenate.sh PA14+SA:A1,A2,A3,B1,B2,B3,C1,C2,C3,D1,D2,D3,E1,E2,E3,F1,F2,F3,G1,G2,G3,H1,H2,H3 ../rawdata/r1r4_welldemultiplex_full_v2/r4/ 4 ../rawdata/condition_demultiplex_full_v2/ &
bash ./concatenate.sh PAO1+BC:A10,A11,A12,B10,B11,B12,C10,C11,C12,D10,D11,D12,E10,E11,E12,F10,F11,F12,G10,G11,G12,H10,H11,H12 ../rawdata/r1r4_welldemultiplex_full_v2/r4/ 4 ../rawdata/condition_demultiplex_full_v2/ &
bash ./concatenate.sh PAO1+SA:A7,A8,A9,B7,B8,B9,C7,C8,C9,D7,D8,D9,E7,E8,E9,F7,F8,F9,G7,G8,G9,H7,H8,H9 ../rawdata/r1r4_welldemultiplex_full_v2/r4/ 4 ../rawdata/condition_demultiplex_full_v2/ &

wait