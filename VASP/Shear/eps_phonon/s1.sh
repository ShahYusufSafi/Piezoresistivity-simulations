source /etc/profile.d/modules.sh 2>/dev/null || true
module load intel-oneapi-mkl

for d in eps_*/; do
    ( cd "$d" || exit
      cp INCAR.relax INCAR
      cp KPOINTS.prim KPOINTS
      echo "=== $d relaxing ==="
      vasp > log.relax 2>&1
      cp CONTCAR POSCAR            # save relaxed cell for the phonon stage
      grep -q "reached required accuracy" OUTCAR && echo "  OK" || echo "  NOT converged"
    )

done
