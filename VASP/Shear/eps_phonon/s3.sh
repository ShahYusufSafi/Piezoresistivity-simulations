export PATH=$HOME/.local/bin:$PATH

for d in eps_*/; do

(
    cd "$d" || exit
    phonopy --fc-symmetry -c POSCAR.relaxed band.conf
    phonopy-bandplot --gnuplot band.yaml > band.dat
    echo "=== outputs: ==="
    ls -l band.yaml band.dat
    echo "=== imaginary (THz)? ==="
    awk 'NF==2 && $2<-0.05{c++} END{print "negatives:",c+0}' band.dat
)
done
