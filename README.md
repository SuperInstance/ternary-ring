# ternary-ring

Ring and field structures for ternary values: Z/3Z arithmetic, polynomial rings, and GF(3^n) field extensions.

## Why This Exists

The integers modulo 3 form a finite field — every nonzero element has a multiplicative inverse. This makes Z/3Z fundamentally different from Z/4Z (which isn't even an integral domain). Finite fields of characteristic 3 (GF(3), GF(9), GF(27), GF(81), ...) are essential in coding theory, cryptography, and combinatorics. This crate provides complete Z/3Z arithmetic, polynomial rings over Z/3Z with GCD and irreducibility testing, and fully functional GF(3^n) field extensions with arithmetic, inversion, order computation, and primitive element generation.

## Core Concepts

- **`Z3`** — An element of Z/3Z (integers modulo 3). Values: `Z3(0)`, `Z3(1)`, `Z3(2)`. Supports `add`, `sub`, `mul`, `neg`, `inv`, `div`.
- **`PolyZ3`** — A polynomial over Z/3Z, stored as coefficients from lowest to highest degree. Supports arithmetic, evaluation, division with remainder, GCD, and irreducibility testing.
- **`GF3n`** — An element of the finite field GF(3^n), represented as a polynomial of degree < n modulo an irreducible polynomial. Full field arithmetic including inversion.

## Quick Start

```toml
# Cargo.toml
[dependencies]
ternary-ring = "0.1"
```

```rust
use ternary_ring::*;

fn main() {
    // Z/3Z arithmetic
    let a = Z3(2);
    let b = Z3(2);
    assert_eq!(a.mul(b), Z3(1)); // 2 * 2 = 4 ≡ 1 (mod 3)
    assert_eq!(a.inv(), Some(Z3(2))); // 2 is its own inverse!

    // Polynomial arithmetic over Z/3Z
    let p = PolyZ3::new(vec![Z3(1), Z3(0), Z3(1)]); // x² + 1
    let q = PolyZ3::new(vec![Z3(1), Z3(1)]);         // x + 1
    let (quotient, remainder) = p.div_rem(&q);
    println!("({:?}) / ({:?}) = ({:?}) rem ({:?})", p, q, quotient, remainder);

    // Irreducibility testing
    assert!(PolyZ3::new(vec![Z3(1), Z3(0), Z3(1)]).is_irreducible()); // x² + 1 over Z/3Z
    assert!(!PolyZ3::new(vec![Z3(1), Z3(2), Z3(1)]).is_irreducible()); // x² + 2x + 1 = (x+1)²

    // GF(9) = GF(3²) field extension
    let field = GF3n::with_degree(2);
    let a = field.element(vec![Z3(1), Z3(1)]); // 1 + x
    let b = field.element(vec![Z3(2), Z3(1)]); // 2 + x
    let c = field.mul(&a, &b);
    let inv = field.inv(&a).unwrap();
    let one = field.mul(&a, &inv);
    println!("a * a⁻¹ = {:?}", one); // should be 1

    // Primitive element (generator of GF(9)*)
    let prim = field.primitive_element();
    println!("Primitive element order: {}", field.element_order(&prim)); // 8 = 3² - 1
}
```

## API Overview

### Z/3Z Arithmetic (`Z3`)
- `Z3::new(v)` / `Z3(val)` — Construct (reduces mod 3)
- `a.add(b)`, `a.sub(b)`, `a.mul(b)`, `a.neg()` — Basic arithmetic
- `a.inv()` — Multiplicative inverse (None for zero)
- `a.div(b)` — Division via inverse
- `Z3::all()` — All three elements [Z3(0), Z3(1), Z3(2)]

### Polynomial Ring (`PolyZ3`)
- `PolyZ3::new(coeffs)` / `zero()` / `constant(c)` / `monomial(c, deg)` — Construction
- `p.degree()` — Degree (None for zero polynomial)
- `p.eval(x)` — Evaluate at a Z/3Z point
- `p.add(&q)`, `p.sub(&q)`, `p.mul(&q)`, `p.scalar_mul(c)` — Arithmetic
- `p.div_rem(&q)` — Polynomial long division → (quotient, remainder)
- `PolyZ3::gcd(&a, &b)` — Euclidean GCD
- `p.is_irreducible()` — Irreducibility test over Z/3Z
- `PolyZ3::irreducibles(degree)` — Enumerate all monic irreducible polynomials of given degree

### Field Extensions (`GF3n`)
- `GF3n::with_degree(n)` — Create GF(3^n) using the first available irreducible
- `GF3n::new(modulus)` — Create with a specific irreducible polynomial
- `field.element(coeffs)` / `field.zero()` / `field.one()` — Construct elements
- `field.add(&a, &b)` / `field.sub(&a, &b)` / `field.mul(&a, &b)` — Field arithmetic
- `field.inv(&a)` / `field.div(&a, &b)` — Inversion and division
- `field.pow(&a, exp)` — Exponentiation by squaring
- `field.element_order(&a)` — Multiplicative order
- `field.primitive_element()` — Find a generator of the multiplicative group
- `field.minimal_polynomial(&a)` — Minimal polynomial over Z/3Z
- `field.field_size()` — 3^n

## How It Works

**Z/3Z arithmetic** is straightforward modular arithmetic. Notably, `2 * 2 = 1 (mod 3)`, so 2 is its own multiplicative inverse.

**Polynomial arithmetic** follows standard rules, with all coefficients reduced mod 3. Division uses polynomial long division with the leading coefficient's inverse.

**Irreducibility testing** first checks for roots in Z/3Z (linear factors), then for higher-degree factors by computing gcd(x^(3^k) − x, p(x)) for k up to n/2 using repeated squaring modulo p(x).

**GF(3^n)** represents elements as polynomials of degree < n. Multiplication multiplies the polynomials and reduces modulo the defining irreducible. Inversion uses the extended Euclidean algorithm to find polynomials s, t such that a(x)s(x) + m(x)t(x) = 1.

## Use Cases

1. **Ternary error-correcting codes** — Build codes over GF(3) and GF(3^n) for ternary channel coding.
2. **Cryptography** — Use GF(3^n) as the basis for discrete-log or pairing-based cryptosystems with characteristic-3 fields.
3. **Combinatorial designs** — Construct difference sets and Latin squares over finite fields of characteristic 3.
4. **Mathematical research** — Enumerate irreducible polynomials, compute minimal polynomials, and explore field structure.

## Ecosystem

- [`ternary-permutation`](https://github.com/user/ternary-permutation) — Permutation groups on ternary vectors
- [`ternary-clustering`](https://github.com/user/ternary-clustering) — Clustering algorithms for ternary data
- [`ternary-automata`](https://github.com/user/ternary-automata) — Cellular automata with ternary states

## License

MIT
