#!/usr/bin/env bash
# Reproduce every result in this repository. Linux or WSL, g++ >= 11, CMake, git. About 30-40 minutes, almost all of it
# building CAPD; the proofs themselves take seconds. Exit status 0 only if every proof run VERIFIES and every
# negative control is REJECTED.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
CAPD_COMMIT=03dc5628203334b214bb7d9fd63788a175521005
CAPD="${CAPD_DIR:-$HERE/_capd}"
if [ ! -f "$CAPD/build/libcapd.a" ]; then
  git clone https://github.com/CAPDGroup/CAPD.git "$CAPD" || exit 1
  git -C "$CAPD" checkout "$CAPD_COMMIT" || exit 1
  cmake -S "$CAPD" -B "$CAPD/build" -DCMAKE_BUILD_TYPE=Release || exit 1
  make -C "$CAPD/build" -j"$(nproc)" || exit 1
fi
echo "CAPD commit: $(git -C "$CAPD" rev-parse HEAD)"
CF=$(sed -n 's/^Cflags: //p' "$CAPD/build/bin/capd.pc"); LB=$(sed -n 's/^Libs: //p' "$CAPD/build/bin/capd.pc")
echo "sha256 of the proof inputs:"; (cd "$HERE" && sha256sum proof/design_chain.txt proof/check_chain.cpp independent/verify_b.cpp)
fail=0
want() {  # want <expected: VERIFIED|REJECTED> <label> <log> <pattern meaning verified>
  if grep -q "$4" "$3"; then got=VERIFIED; else got=REJECTED; fi
  if [ "$got" = "$1" ]; then echo "  OK    $2: $got"; else echo "  WRONG $2: $got (expected $1)"; fail=1; fi
}

cd "$HERE/proof" || exit 1
g++ $CF check_chain.cpp $LB -o check_chain || exit 1
g++ $CF -DUEDA_NEGATIVE_CONTROL_BUILD check_chain.cpp $LB -o controls/check_chain_negctl || exit 1
g++ $CF controls/cmp_test.cpp $LB -o controls/cmp_test || exit 1
g++ $CF "$CAPD/capdDynSys/examples/RosslerChaoticDynamics/RosslerChaoticDynamics.cpp" $LB -o controls/rossler_control || exit 1

echo "== control 1: CAPD interval comparisons mean 'certainly'"
./controls/cmp_test > controls/cmp_test.log 2>&1
want VERIFIED "comparison semantics" controls/cmp_test.log "^COMPARISON SEMANTICS OK"
echo "== control 2: CAPD's own proof of the Rossler horseshoe (Zgliczynski 1997); every line must read true"
./controls/rossler_control | tee controls/rossler_control.log
grep -q false controls/rossler_control.log && { echo "  WRONG Rossler control"; fail=1; } || echo "  OK    Rossler control"

echo "== checker A on proof/design_chain.txt (Taylor orders 12, 10, 14)"
{ time ./check_chain design_chain.txt 64 40; } > run_A.log 2>&1
UEDA_ORDER=10 ./check_chain design_chain.txt 64 40 > run_A_order10.log 2>&1
UEDA_ORDER=14 ./check_chain design_chain.txt 64 40 > run_A_order14.log 2>&1
want VERIFIED "checker A, order 12" run_A.log "^ALL COVERING RELATIONS AND DISJOINTNESS VERIFIED"
want VERIFIED "checker A, order 10" run_A_order10.log "^ALL COVERING RELATIONS AND DISJOINTNESS VERIFIED"
want VERIFIED "checker A, order 14" run_A_order14.log "^ALL COVERING RELATIONS AND DISJOINTNESS VERIFIED"
grep -E "^(disjoint|set|transition)" run_A.log | sed 's/^/    /'

echo "== negative controls, checker A (each must be REJECTED)"
for d in neg_orient neg_narrowC3 neg_thinC1; do
  ./check_chain controls/$d.txt 64 40 > controls/run_A_$d.log 2>&1
  want REJECTED "checker A, $d" controls/run_A_$d.log "^ALL COVERING RELATIONS AND DISJOINTNESS VERIFIED"
done
for b in 7.0 7.4; do
  UEDA_NEGCTL_B=$b ./controls/check_chain_negctl design_chain.txt 64 40 > controls/run_A_B$b.log 2>&1
  want REJECTED "checker A, B = $b" controls/run_A_B$b.log "^ALL COVERING RELATIONS AND DISJOINTNESS VERIFIED"
done

cd "$HERE/independent" || exit 1
g++ $CF verify_b.cpp $LB -o verify_b || exit 1
echo "== checker B on proof/design_chain.txt"
{ time ./verify_b ../proof/design_chain.txt; } > run_B.log 2>&1
want VERIFIED "checker B" run_B.log "^ALL VERIFIED"
echo "== negative controls, checker B (each must be REJECTED)"
./verify_b ../proof/design_chain.txt --flip-last > run_B_flip.log 2>&1
want REJECTED "checker B, last orientation flipped" run_B_flip.log "^ALL VERIFIED"
./verify_b ../proof/design_chain.txt --B 7.4 > run_B_B7.4.log 2>&1
want REJECTED "checker B, B = 7.4" run_B_B7.4.log "^ALL VERIFIED"

echo
if [ $fail -eq 0 ]; then echo "REPRODUCE: ALL RESULTS AND CONTROLS AS STATED"; else echo "REPRODUCE: MISMATCH (see above)"; fi
exit $fail
