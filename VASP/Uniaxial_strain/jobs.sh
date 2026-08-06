#!/bin/bash
set -euo pipefail

USER="ws2505"
HOST="stud2.mpi.univie.ac.at"
REMOTE="~/Ysafi_strian/uniaxial"

read -s -p "Password: " PASS; echo

SSH="sshpass -p ${PASS} ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null ${USER}@${HOST}"
SCP="sshpass -p ${PASS} scp -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null"


cat > jobscript << 'EOF'
source /etc/profile.d/modules.sh 2>/dev/null || true
module load intel-oneapi-mkl

# Stage A - ionic relaxation (shear only; skipped if INCAR.relax absent)
if [ -f INCAR.relax ]; then
    cp INCAR.relax INCAR
    cp KPOINTS.mesh KPOINTS
    vasp > log.relax
    cp CONTCAR POSCAR
fi

# Stage B - SCF on uniform mesh -> CHGCAR  (STRESS TENSOR LIVES HERE)
cp INCAR.scf INCAR      || { echo "Missing INCAR.scf";    exit 1; }
cp KPOINTS.mesh KPOINTS || { echo "Missing KPOINTS.mesh"; exit 1; }
vasp > log.scf
grep -q "reached required accuracy" log.scf || echo "WARNING: SCF not converged"
cp OUTCAR OUTCAR.scf ; cp vasprun.xml vasprun.scf.xml     # <-- preserve before Stage C

# Stage C - non-SCF bands on the k-path, reads CHGCAR
cp INCAR.band INCAR     || { echo "Missing INCAR.band";   exit 1; }
cp KPOINTS.line KPOINTS || { echo "Missing KPOINTS.line"; exit 1; }
vasp > log.band
cp OUTCAR OUTCAR.band ; cp vasprun.xml vasprun.band.xml
EOF



# Phase 1: upload and run each strain
for d in eps_*/; do
    cp jobscript ${d}
    
    name="${d%/}"
    echo "=== $name ==="
    $SSH "rm -rf ${REMOTE}/${name} && mkdir -p ${REMOTE}/${name}"   # 2. wipe+recreate
    #$SCP ${d}{INCAR.relax,INCAR.scf,INCAR.band,KPOINTS.mesh,KPOINTS.line,POSCAR,POTCAR,jobscript} \
    #    ${USER}@${HOST}:${REMOTE}/${name}/ 2>/dev/null || true

    $SCP ${d}{INCAR.scf,INCAR.band,KPOINTS.mesh,KPOINTS.line,POSCAR,POTCAR,jobscript} \
        ${USER}@${HOST}:${REMOTE}/${name}/ 2>/dev/null || true       # 3. upload (jobscript incl.)
    echo "files copied"
    $SSH "cd ${REMOTE}/${name} && bash --login jobscript"
    echo "  done"
done

# Phase 2: retrieve outputs
for d in eps_*/; do
    name="${d%/}"
    echo "=== retrieving $name ==="
    for f in OUTCAR.scf OUTCAR.band vasprun.scf.xml vasprun.band.xml OSZICAR CONTCAR log.relax log.scf log.band; do
        $SCP ${USER}@${HOST}:${REMOTE}/${name}/${f} ${d} \
            2>/dev/null || echo "  missing: ${f}"
    done
done

echo "=== all done ==="