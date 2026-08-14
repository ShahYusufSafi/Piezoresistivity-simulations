#!/bin/bash
set -euo pipefail

USER="ws2505"
HOST="stud2.mpi.univie.ac.at"
REMOTE="/home/ws2505/Ysafi_runs/Strain_Sweeps/shear"

read -s -p "Password: " PASS; echo
SSH="sshpass -p ${PASS} ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null ${USER}@${HOST}"
SCP="sshpass -p ${PASS} scp -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null"

STRAIN="eps_+0.0000"
SRC="${STRAIN}/mass_val"                    # local: inputs here, outputs land here
RUNDIR="${REMOTE}/${STRAIN}/mass_val"       # remote: isolated mass subdir

cat > ${SRC}/jobscript << 'EOF'
source /etc/profile.d/modules.sh 2>/dev/null || true
module load intel-oneapi-mkl
cp INCAR.band INCAR
cp KPOINTS.grid KPOINTS
vasp > log.mass 2>&1
grep -q "General timing" OUTCAR && echo "VASP OK" || echo "VASP FAILED"
EOF

# isolated mass dir on the server; upload everything INCLUDING CHGCAR
$SSH "mkdir -p ${RUNDIR}"
$SCP ${SRC}/{POSCAR,POTCAR,INCAR.band,KPOINTS.grid,CHGCAR,jobscript} \
     ${USER}@${HOST}:${RUNDIR}/

# run
$SSH "cd ${RUNDIR} && bash --login jobscript"

# retrieve into the same local mass_val
for f in OUTCAR OSZICAR vasprun.xml EIGENVAL log.mass; do
    $SCP ${USER}@${HOST}:${RUNDIR}/${f} ${SRC}/ 2>/dev/null || echo "missing: ${f}"
done
echo "=== done, outputs in ${SRC}/ ==="