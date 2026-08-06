#!/bin/bash
set -euo pipefail

USER="ws2505"
HOST="stud2.mpi.univie.ac.at"
REMOTE="~/Ysafi_strian/shear"

read -s -p "Password: " PASS; echo

SSH="sshpass -p ${PASS} ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null ${USER}@${HOST}"
SCP="sshpass -p ${PASS} scp -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null"

cat > jobscript << 'EOF'
source /etc/profile.d/modules.sh 2>/dev/null || true
module load intel-oneapi-mkl

# Stage A - ionic relaxation (Kleinman displacement under [111] shear)
cp INCAR.relax INCAR; cp KPOINTS.mesh KPOINTS
vasp > log.relax
cp CONTCAR POSCAR

# Stage B - SCF on uniform mesh -> writes CHGCAR
cp INCAR.scf INCAR; cp KPOINTS.mesh KPOINTS
vasp > log.scf

# Stage C - non-SCF bands on the k-path, reads CHGCAR
cp INCAR.band INCAR; cp KPOINTS.line KPOINTS
vasp > log.band
EOF

# Phase 1: upload and run each strain
for d in eps_*/; do
    cp jobscript ${d}
    name="${d%/}"
    echo "=== $name ==="
    $SSH "rm -rf ${REMOTE}/${name} && mkdir -p ${REMOTE}/${name}"
    $SCP ${d}{INCAR.relax,INCAR.scf,INCAR.band,KPOINTS.mesh,KPOINTS.line,POSCAR,POTCAR,jobscript} \
        ${USER}@${HOST}:${REMOTE}/${name}/
    echo "files copied"
    $SSH "cd ${REMOTE}/${name} && bash --login jobscript"
    echo "  done"
done

# Phase 2: retrieve outputs
for d in eps_*/; do
    name="${d%/}"
    echo "=== retrieving $name ==="
    for f in OUTCAR OSZICAR vasprun.xml CONTCAR log.relax log.scf log.band; do
        $SCP ${USER}@${HOST}:${REMOTE}/${name}/${f} ${d} \
            2>/dev/null || echo "  missing: ${f}"
    done
done

echo "=== all done ==="