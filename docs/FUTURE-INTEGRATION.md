# Future Integration: ternary-ring

## Current State
Implements ring and field structures for ternary arithmetic: `Z3` (integers modulo 3 with full arithmetic), `PolyZ3` (polynomial rings over Z/3Z with Euclidean division), GF(3^n) field extensions for n=2,3,4, minimal polynomial computation, irreducibility testing, and polynomial evaluation.

## Integration Opportunities

### With ternary-matrix
Matrix operations over GF(3) instead of over the reals. `TernaryMatrix::gf3_inverse()` already computes inverses over GF(3) — this crate provides the mathematical foundation. `PolyZ3` enables characteristic polynomial computation over GF(3), giving exact eigenvalues (no floating-point approximation). The minimal polynomial of a matrix over GF(3) determines its rational canonical form.

### With ternary-codes
Error-correcting codes over GF(3). `PolyZ3` with `irreducible()` testing constructs generator polynomials for ternary BCH codes. `GF(3^n)` field extensions provide the algebraic structure for Reed-Solomon codes over ternary alphabets. `poly_eval()` encodes messages; the Euclidean algorithm (`poly_div_rem()`) decodes them.

### With ternary-cell / construct-core
Room state algebra. Room states are vectors over Z/3Z. State transitions are linear maps over GF(3). `Z3::add()` and `Z3::mul()` define the arithmetic for cell state combination. Conservation laws are linear equations over Z/3Z — solvable exactly via GF(3) linear algebra. No floating-point rounding, ever.

## Potential in Mature Systems
In PLATO, `Z3` and `PolyZ3` are the fundamental algebraic types. Every computation in the system reduces to operations over Z/3Z. The `GF(3^n)` extensions enable sophisticated algebraic constructions at Layer 1: error-correcting codes for reliable communication, polynomial commitment schemes for construct integrity verification, and finite field arithmetic for cryptographic operations. At Layer 0, `Z3` arithmetic compiles to a handful of lookup tables — addition and multiplication tables are 9 bytes each.

## Cross-Pollination Ideas
**Music × Ring:** Ternary pitch class sets form a module over Z/3Z. Transposition is addition; inversion is negation. `PolyZ3` evaluated at specific points generates chord voicings. The irreducible polynomials over GF(3) correspond to "indecomposable" harmonic structures — chords that cannot be factored into simpler progressions. This is algebraic music theory. Connects to `flux-algebra-rs`.

**Cryptography × Ring:** GF(3^n) is the foundation for ternary lattice-based cryptography. `PolyZ3` polynomial multiplication (using NTT over GF(3)) enables fast lattice operations. `irreducible()` testing constructs good lattice bases. This is post-quantum cryptography on ternary hardware.

**Physics × Ring:** Gauge theories over GF(3) model ternary physical systems. The `PolyZ3` ring provides the algebra for constructing gauge-invariant observables. Conservation laws in construct-core are Noether currents over this ring.

## Dependencies for Next Steps
- NTT (Number Theoretic Transform) over GF(3) for fast polynomial multiplication
- Larger field extensions (GF(3^8), GF(3^16)) for cryptographic applications
- Integration with `ternary-matrix` for GF(3) linear algebra
- Connection to `ternary-codes` for error-correcting code construction
