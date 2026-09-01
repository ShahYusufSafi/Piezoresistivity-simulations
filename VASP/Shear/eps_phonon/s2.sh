source /etc/profile.d/modules.sh 2>/dev/null || true
module load intel-oneapi-mkl
export PATH=$HOME/.local/bin:$PATH

echo "=== phonopy version ==="
phonopy --version

for d in eps_*/; do
    (
    cd "$d" || exit
    cp INCAR.force INCAR
    cp KPOINTS.sc KPOINTS

    phonopy -d --dim 2 2 2 -c POSCAR.relaxed

    NDISP=$(ls POSCAR-0* 2>/dev/null | wc -l)
    echo "=== displacements found: $NDISP ==="
    [ "$NDISP" -eq 0 ] && { echo "still zero, stopping"; exit 1; }
    ls POSCAR-0* # debugging 

    echo "=== forces ==="
    for f in POSCAR-0*; do
        n=${f#POSCAR-}
        cp "$f" POSCAR
        vasp > log.force-$n 2>&1
        cp vasprun.xml vasprun-$n.xml
        echo "  done: $n"
    done

    echo "=== FORCE_SETS ==="
    phonopy -f vasprun-*.xml
    ls -l FORCE_SETS && echo SUCCESS for eps $d || echo FAILED

    )
done
