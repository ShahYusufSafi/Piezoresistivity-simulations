source /etc/profile.d/modules.sh 2>/dev/null || true
module load intel-oneapi-mkl


for d in eps_*/; do
    ( cd "$d" || exit

      #rm -f OUTCAR CONTCAR WAVECAR CHG CHGCAR vasprun.xml OSZICAR log.* POSCAR.relaxed       #IBZKPT EIGENVAL DOSCAR PCDAT XDATCAR REPORT

      cp INCAR.relax INCAR
      cp KPOINTS.prim KPOINTS
      echo "=== $d relaxing ==="
      vasp > log.relax 2>&1
      cp CONTCAR POSCAR.relaxed            # save relaxed cell for the phonon stage
      grep -q "reached required accuracy" OUTCAR && echo "  OK" || echo "  NOT converged"
    )
done
