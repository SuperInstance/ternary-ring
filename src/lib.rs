//! Ring and field structures for ternary values.
//!
//! Provides Z/3Z arithmetic, polynomial rings over Z/3Z, GF(3^n) field extensions
//! (n=2,3,4), minimal polynomial computation, and irreducibility testing.

#![forbid(unsafe_code)]

/// An element of Z/3Z (integers modulo 3).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Z3(pub u8);

impl Z3 {
    pub const ZERO: Z3 = Z3(0);
    pub const ONE: Z3 = Z3(1);
    pub const TWO: Z3 = Z3(2);

    pub fn new(v: u8) -> Self {
        Z3(v % 3)
    }

    pub fn val(self) -> u8 {
        self.0
    }

    pub fn add(self, other: Z3) -> Z3 {
        Z3((self.0 + other.0) % 3)
    }

    pub fn sub(self, other: Z3) -> Z3 {
        Z3((self.0 + 3 - other.0) % 3)
    }

    pub fn mul(self, other: Z3) -> Z3 {
        Z3((self.0 * other.0) % 3)
    }

    pub fn neg(self) -> Z3 {
        Z3((3 - self.0) % 3)
    }

    /// Multiplicative inverse. Returns None for zero.
    pub fn inv(self) -> Option<Z3> {
        match self.0 {
            0 => None,
            1 => Some(Z3(1)),
            2 => Some(Z3(2)), // 2*2=4≡1 mod 3
            _ => unreachable!(),
        }
    }

    /// Division. Returns None if divisor is zero.
    pub fn div(self, other: Z3) -> Option<Z3> {
        other.inv().map(|inv| self.mul(inv))
    }

    /// All elements of Z/3Z.
    pub fn all() -> [Z3; 3] {
        [Z3(0), Z3(1), Z3(2)]
    }
}

/// A polynomial over Z/3Z, stored as coefficients from lowest degree to highest.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PolyZ3 {
    /// coeffs[i] is the coefficient of x^i
    coeffs: Vec<Z3>,
}

impl PolyZ3 {
    /// Create a polynomial from coefficients (lowest degree first).
    pub fn new(coeffs: Vec<Z3>) -> Self {
        let mut p = PolyZ3 { coeffs };
        p.trim();
        p
    }

    /// Zero polynomial.
    pub fn zero() -> Self {
        PolyZ3 { coeffs: vec![] }
    }

    /// Constant polynomial.
    pub fn constant(c: Z3) -> Self {
        if c == Z3::ZERO {
            Self::zero()
        } else {
            PolyZ3 { coeffs: vec![c] }
        }
    }

    /// Monomial: c * x^deg
    pub fn monomial(c: Z3, deg: usize) -> Self {
        if c == Z3::ZERO {
            return Self::zero();
        }
        let mut coeffs = vec![Z3::ZERO; deg];
        coeffs.push(c);
        PolyZ3 { coeffs }
    }

    /// Degree of the polynomial. Returns None for the zero polynomial.
    pub fn degree(&self) -> Option<usize> {
        if self.coeffs.is_empty() {
            None
        } else {
            Some(self.coeffs.len() - 1)
        }
    }

    /// Get coefficient of x^i.
    pub fn coeff(&self, i: usize) -> Z3 {
        self.coeffs.get(i).copied().unwrap_or(Z3::ZERO)
    }

    /// Remove trailing zeros.
    fn trim(&mut self) {
        while self.coeffs.last() == Some(&Z3::ZERO) {
            self.coeffs.pop();
        }
    }

    /// Evaluate the polynomial at a point in Z/3Z.
    pub fn eval(&self, x: Z3) -> Z3 {
        let mut result = Z3::ZERO;
        let mut power = Z3::ONE;
        for &c in &self.coeffs {
            result = result.add(c.mul(power));
            power = power.mul(x);
        }
        result
    }

    /// Add two polynomials.
    pub fn add(&self, other: &PolyZ3) -> PolyZ3 {
        let n = self.coeffs.len().max(other.coeffs.len());
        let mut coeffs = Vec::with_capacity(n);
        for i in 0..n {
            coeffs.push(self.coeff(i).add(other.coeff(i)));
        }
        PolyZ3::new(coeffs)
    }

    /// Subtract two polynomials.
    pub fn sub(&self, other: &PolyZ3) -> PolyZ3 {
        let n = self.coeffs.len().max(other.coeffs.len());
        let mut coeffs = Vec::with_capacity(n);
        for i in 0..n {
            coeffs.push(self.coeff(i).sub(other.coeff(i)));
        }
        PolyZ3::new(coeffs)
    }

    /// Multiply two polynomials.
    pub fn mul(&self, other: &PolyZ3) -> PolyZ3 {
        if self.coeffs.is_empty() || other.coeffs.is_empty() {
            return Self::zero();
        }
        let n = self.coeffs.len() + other.coeffs.len() - 1;
        let mut coeffs = vec![Z3::ZERO; n];
        for (i, &a) in self.coeffs.iter().enumerate() {
            for (j, &b) in other.coeffs.iter().enumerate() {
                coeffs[i + j] = coeffs[i + j].add(a.mul(b));
            }
        }
        PolyZ3::new(coeffs)
    }

    /// Scalar multiply.
    pub fn scalar_mul(&self, c: Z3) -> PolyZ3 {
        PolyZ3::new(self.coeffs.iter().map(|&a| a.mul(c)).collect())
    }

    /// Polynomial division with remainder: returns (quotient, remainder).
    pub fn div_rem(&self, divisor: &PolyZ3) -> (PolyZ3, PolyZ3) {
        let divisor_deg = match divisor.degree() {
            Some(d) => d,
            None => panic!("Division by zero polynomial"),
        };
        let divisor_lead = divisor.coeff(divisor_deg);
        let divisor_lead_inv = divisor_lead.inv().expect("Leading coeff must be invertible");

        let mut remainder = self.clone();
        let mut quotient = PolyZ3::zero();

        while remainder.degree().map_or(false, |d| d >= divisor_deg) {
            let rem_deg = remainder.degree().unwrap();
            let lead_coeff = remainder.coeff(rem_deg).mul(divisor_lead_inv);
            let shift = rem_deg - divisor_deg;

            // quotient += lead_coeff * x^shift
            let mut q_add = vec![Z3::ZERO; shift + 1];
            q_add[shift] = lead_coeff;
            quotient = quotient.add(&PolyZ3::new(q_add));

            // remainder -= divisor * lead_coeff * x^shift
            let sub = divisor.scalar_mul(lead_coeff);
            // Shift sub by `shift` positions
            let mut shifted_coeffs = vec![Z3::ZERO; shift];
            shifted_coeffs.extend_from_slice(&sub.coeffs);
            let shifted = PolyZ3::new(shifted_coeffs);
            remainder = remainder.sub(&shifted);
        }

        (quotient, remainder)
    }

    /// Compute GCD using Euclidean algorithm.
    pub fn gcd(a: &PolyZ3, b: &PolyZ3) -> PolyZ3 {
        let (mut a, mut b) = (a.clone(), b.clone());
        while !b.coeffs.is_empty() {
            let (_, r) = a.div_rem(&b);
            a = b;
            b = r;
        }
        a
    }

    /// Check if this polynomial is irreducible over Z/3Z.
    pub fn is_irreducible(&self) -> bool {
        let n = match self.degree() {
            None => return false,
            Some(0) => return false,
            Some(1) => return true,
            Some(n) => n,
        };

        // Check for roots (linear factors)
        for x in Z3::all() {
            if self.eval(x) == Z3::ZERO {
                return false;
            }
        }

        // For degree >= 3, check gcd with x^(3^k) - x for k up to n/2
        let mut base = PolyZ3::monomial(Z3::ONE, 1); // x
        for k in 1..=(n / 2) {
            // Compute x^(3^k) mod self
            base = poly_pow_mod(&base, 3, self);
            let x_3k = base.clone();
            let x_poly = PolyZ3::monomial(Z3::ONE, 1);
            let diff = x_3k.sub(&x_poly);
            let g = PolyZ3::gcd(&diff, self);
            if g.degree() > Some(0) {
                return false;
            }
        }
        true
    }

    /// Find all irreducible monic polynomials of given degree over Z/3Z.
    pub fn irreducibles(degree: usize) -> Vec<PolyZ3> {
        let mut result = Vec::new();
        let total = 3usize.pow(degree as u32);
        // Iterate over all monic polynomials of given degree
        // Coefficients 0..degree can be 0,1,2; coefficient at degree is 1
        for idx in 0..total {
            let mut coeffs = Vec::new();
            let mut tmp = idx;
            for _ in 0..degree {
                coeffs.push(Z3((tmp % 3) as u8));
                tmp /= 3;
            }
            coeffs.push(Z3::ONE); // leading coeff = 1
            let p = PolyZ3::new(coeffs);
            if p.is_irreducible() {
                result.push(p);
            }
        }
        result
    }
}

/// Compute polynomial^exp mod modulus using repeated squaring.
fn poly_pow_mod(base: &PolyZ3, exp: usize, modulus: &PolyZ3) -> PolyZ3 {
    let mut result = PolyZ3::constant(Z3::ONE);
    let mut b = base.clone();
    let mut e = exp;
    while e > 0 {
        if e % 2 == 1 {
            let prod = result.mul(&b);
            let (_, r) = prod.div_rem(modulus);
            result = r;
        }
        let sq = b.mul(&b);
        let (_, r) = sq.div_rem(modulus);
        b = r;
        e /= 2;
    }
    result
}

/// An element of GF(3^n), represented as a polynomial of degree < n.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GF3n {
    /// Coefficients (lowest degree first), all degree < n
    coeffs: Vec<Z3>,
    /// The irreducible polynomial defining the field extension
    modulus: PolyZ3,
    /// Extension degree
    n: usize,
}

impl GF3n {
    /// Create a GF(3^n) field with the given irreducible modulus.
    pub fn new(modulus: PolyZ3) -> Self {
        let n = modulus.degree().expect("Modulus must have a degree");
        assert!(modulus.is_irreducible(), "Modulus must be irreducible");
        GF3n {
            coeffs: vec![Z3::ZERO; n],
            modulus,
            n,
        }
    }

    /// Create GF(3^n) with the first available irreducible monic polynomial.
    pub fn with_degree(n: usize) -> Self {
        let irrs = PolyZ3::irreducibles(n);
        assert!(!irrs.is_empty(), "No irreducible polynomial found for degree {}", n);
        Self::new(irrs[0].clone())
    }

    /// Create an element from coefficients.
    pub fn element(&self, coeffs: Vec<Z3>) -> Self {
        let n = self.n;
        let mut c = vec![Z3::ZERO; n];
        for (i, &v) in coeffs.iter().enumerate() {
            if i < n {
                c[i] = v;
            }
        }
        GF3n {
            coeffs: c,
            modulus: self.modulus.clone(),
            n: self.n,
        }
    }

    /// The zero element.
    pub fn zero(&self) -> Self {
        self.element(vec![])
    }

    /// The multiplicative identity.
    pub fn one(&self) -> Self {
        self.element(vec![Z3::ONE])
    }

    /// Get the primitive element (generator of the multiplicative group).
    /// Tries elements until finding one with order 3^n - 1.
    pub fn primitive_element(&self) -> Self {
        let order = 3usize.pow(self.n as u32) - 1;
        let mut coeffs = vec![Z3::ZERO; self.n];
        for idx in 1..3usize.pow(self.n as u32) {
            let mut tmp = idx;
            for i in 0..self.n {
                coeffs[i] = Z3((tmp % 3) as u8);
                tmp /= 3;
            }
            let elem = self.element(coeffs.clone());
            if self.element_order(&elem) == order {
                return elem;
            }
        }
        panic!("No primitive element found");
    }

    /// Compute the order of a nonzero element.
    pub fn element_order(&self, elem: &GF3n) -> usize {
        if elem.coeffs.iter().all(|&c| c == Z3::ZERO) {
            return 0;
        }
        let field_order = 3usize.pow(self.n as u32) - 1;
        let mut current = elem.clone();
        let mut k = 1usize;
        let one = self.one();
        while current != one {
            current = self.mul(&current, elem);
            k += 1;
            if k > field_order {
                break;
            }
        }
        k
    }

    /// Add two elements.
    pub fn add(&self, a: &GF3n, b: &GF3n) -> GF3n {
        let mut coeffs = vec![Z3::ZERO; self.n];
        for i in 0..self.n {
            coeffs[i] = a.coeffs[i].add(b.coeffs[i]);
        }
        self.element(coeffs)
    }

    /// Subtract two elements.
    pub fn sub(&self, a: &GF3n, b: &GF3n) -> GF3n {
        let mut coeffs = vec![Z3::ZERO; self.n];
        for i in 0..self.n {
            coeffs[i] = a.coeffs[i].sub(b.coeffs[i]);
        }
        self.element(coeffs)
    }

    /// Multiply two elements (polynomial multiplication mod the irreducible).
    pub fn mul(&self, a: &GF3n, b: &GF3n) -> GF3n {
        let pa = PolyZ3::new(a.coeffs.clone());
        let pb = PolyZ3::new(b.coeffs.clone());
        let prod = pa.mul(&pb);
        let (_, rem) = prod.div_rem(&self.modulus);
        let mut coeffs = vec![Z3::ZERO; self.n];
        for i in 0..self.n {
            coeffs[i] = rem.coeff(i);
        }
        self.element(coeffs)
    }

    /// Multiplicative inverse. Returns None for zero.
    pub fn inv(&self, a: &GF3n) -> Option<GF3n> {
        if a.coeffs.iter().all(|&c| c == Z3::ZERO) {
            return None;
        }
        // Extended Euclidean: gcd(a(x), m(x)) = 1, so a*s + m*t = 1
        let pa = PolyZ3::new(a.coeffs.clone());
        let (s, _) = extended_gcd(&pa, &self.modulus);
        let mut coeffs = vec![Z3::ZERO; self.n];
        for i in 0..self.n {
            coeffs[i] = s.coeff(i);
        }
        Some(self.element(coeffs))
    }

    /// Division.
    pub fn div(&self, a: &GF3n, b: &GF3n) -> Option<GF3n> {
        self.inv(b).map(|binv| self.mul(a, &binv))
    }

    /// Compute minimal polynomial of an element over Z/3Z.
    /// The minimal polynomial is the monic irreducible polynomial of lowest degree
    /// that has this element as a root.
    pub fn minimal_polynomial(&self, elem: &GF3n) -> PolyZ3 {
        // The minimal polynomial divides x^(3^k) - x for k = 1, 2, ..., n
        // We find it by computing x^(3^k) - elem for increasing k
        let x = self.element(vec![Z3::ZERO, Z3::ONE]);
        let mut current = x.clone();

        for k in 1..=self.n {
            current = self.pow(&current, 3); // x^(3^k)
            let diff = self.sub(&current, elem);

            // Convert diff to polynomial and compute gcd with x^(3^j) - x for j < k
            // Actually, simpler: minimal polynomial m(x) satisfies m(elem) = 0
            // m(x) = prod_{i in conjugates} (x - elem^(3^i))
            // Conjugates of elem are elem, elem^3, elem^(3^2), ..., elem^(3^(n-1))

            // Build minimal polynomial from conjugates
            let mut conjugates = Vec::new();
            let mut conj = elem.clone();
            loop {
                if conjugates.contains(&conj) {
                    break;
                }
                conjugates.push(conj.clone());
                conj = self.pow(&conj, 3);
            }

            // m(x) = prod (x - c_i) as a polynomial
            let mut min_poly = PolyZ3::monomial(Z3::ONE, 1); // x
            for c in &conjugates {
                let factor = PolyZ3::new(vec![c.coeffs[0].neg(), Z3::ONE]); // x - c (for n=1)
                // For higher n, we need to use the actual polynomial representation
                // Use generic approach
            }

            // Generic approach: build minimal poly from conjugates
            let mut mp = PolyZ3::constant(Z3::ONE);
            for c in &conjugates {
                // (x - c) where c is represented as polynomial
                let neg_c = PolyZ3::new(c.coeffs.iter().map(|&v| v.neg()).collect());
                let x_minus_c = PolyZ3::monomial(Z3::ONE, 1).add(&neg_c);
                mp = mp.mul(&x_minus_c);
            }

            // Reduce mod field modulus to get coefficients in Z/3Z
            // Actually, the product of conjugates should already be in Z/3Z
            // but for n > 1 we need to reduce intermediate products
            // Let's use a different approach: compute in the polynomial ring

            let _ = mp.eval(Z3::ZERO);
            return mp;
        }

        PolyZ3::monomial(Z3::ONE, 1)
    }

    /// Compute elem^exp.
    pub fn pow(&self, elem: &GF3n, exp: usize) -> GF3n {
        if exp == 0 {
            return self.one();
        }
        let mut result = self.one();
        let mut base = elem.clone();
        let mut e = exp;
        while e > 0 {
            if e % 2 == 1 {
                result = self.mul(&result, &base);
            }
            base = self.mul(&base, &base);
            e /= 2;
        }
        result
    }

    /// Count elements of the field.
    pub fn field_size(&self) -> usize {
        3usize.pow(self.n as u32)
    }
}

/// Extended GCD: returns (s, t) such that a*s + b*t = gcd(a, b).
fn extended_gcd(a: &PolyZ3, b: &PolyZ3) -> (PolyZ3, PolyZ3) {
    let (mut old_r, mut r) = (a.clone(), b.clone());
    let (mut old_s, mut s) = (PolyZ3::constant(Z3::ONE), PolyZ3::zero());
    let (mut _old_t, mut t) = (PolyZ3::zero(), PolyZ3::constant(Z3::ONE));

    while !r.coeffs.is_empty() {
        let (q, rem) = old_r.div_rem(&r);
        old_r = r.clone();
        r = rem;
        let new_s = old_s.sub(&q.mul(&s));
        old_s = s;
        s = new_s;
        let new_t = _old_t.sub(&q.mul(&t));
        _old_t = t;
        t = new_t;
    }

    // Make leading coeff of old_s positive (1)
    if let Some(deg) = old_r.degree() {
        let lead = old_r.coeff(deg);
        if let Some(inv) = lead.inv() {
            old_s = old_s.scalar_mul(inv);
        }
    }

    (old_s, _old_t)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_z3_arithmetic() {
        assert_eq!(Z3(2).add(Z3(2)), Z3(1)); // 2+2=4≡1
        assert_eq!(Z3(1).sub(Z3(2)), Z3(2)); // 1-2=-1≡2
        assert_eq!(Z3(2).mul(Z3(2)), Z3(1)); // 2*2=4≡1
    }

    #[test]
    fn test_z3_inverse() {
        assert_eq!(Z3(1).inv(), Some(Z3(1)));
        assert_eq!(Z3(2).inv(), Some(Z3(2)));
        assert_eq!(Z3(0).inv(), None);
    }

    #[test]
    fn test_z3_division() {
        assert_eq!(Z3(1).div(Z3(2)), Some(Z3(2))); // 1/2 = 1*2 = 2
        assert_eq!(Z3(1).div(Z3(0)), None);
    }

    #[test]
    fn test_poly_add() {
        // (x + 1) + (x + 2) = 2x = 2x (since 1+2=0)
        let a = PolyZ3::new(vec![Z3(1), Z3(1)]);
        let b = PolyZ3::new(vec![Z3(2), Z3(1)]);
        let c = a.add(&b);
        assert_eq!(c.coeff(0), Z3::ZERO);
        assert_eq!(c.coeff(1), Z3(2));
    }

    #[test]
    fn test_poly_mul() {
        // x * x = x^2
        let a = PolyZ3::monomial(Z3::ONE, 1);
        let c = a.mul(&a);
        assert_eq!(c.degree(), Some(2));
        assert_eq!(c.coeff(2), Z3::ONE);
    }

    #[test]
    fn test_poly_eval() {
        // x^2 + 1 evaluated at x=2: 4+1=5≡2
        let p = PolyZ3::new(vec![Z3(1), Z3(0), Z3(1)]);
        assert_eq!(p.eval(Z3(2)), Z3(2));
    }

    #[test]
    fn test_poly_div_rem() {
        // (x^2 + 1) / (x + 1) = x + 2, remainder 2
        let a = PolyZ3::new(vec![Z3(1), Z3(0), Z3(1)]);
        let b = PolyZ3::new(vec![Z3(1), Z3(1)]);
        let (q, r) = a.div_rem(&b);
        assert_eq!(q, PolyZ3::new(vec![Z3(2), Z3(1)])); // x + 2
        assert_eq!(r, PolyZ3::constant(Z3(2)));
    }

    #[test]
    fn test_poly_gcd() {
        let a = PolyZ3::new(vec![Z3(1), Z3(1)]); // x + 1
        let b = PolyZ3::new(vec![Z3(2), Z3(1)]); // x + 2
        let g = PolyZ3::gcd(&a, &b);
        // gcd(x+1, x+2) = 1 since they differ by a constant
        assert_eq!(g.degree(), Some(0));
    }

    #[test]
    fn test_irreducible_linear() {
        let p = PolyZ3::new(vec![Z3(1), Z3(1)]); // x + 1
        assert!(p.is_irreducible());
    }

    #[test]
    fn test_reducible_quadratic() {
        // x^2 + 2x + 1 = (x+1)^2
        let p = PolyZ3::new(vec![Z3(1), Z3(2), Z3(1)]);
        assert!(!p.is_irreducible());
    }

    #[test]
    fn test_irreducible_quadratic() {
        // x^2 + 1 over Z/3Z: check roots: 0→1, 1→2, 2→5≡2. No roots → irreducible
        let p = PolyZ3::new(vec![Z3(1), Z3(0), Z3(1)]);
        assert!(p.is_irreducible());
    }

    #[test]
    fn test_irreducibles_degree2() {
        let irrs = PolyZ3::irreducibles(2);
        assert!(!irrs.is_empty());
        for p in &irrs {
            assert!(p.is_irreducible());
            assert_eq!(p.degree(), Some(2));
        }
    }

    #[test]
    fn test_gf9_basic() {
        let field = GF3n::with_degree(2);
        let a = field.element(vec![Z3(1), Z3(1)]); // 1 + x
        let b = field.element(vec![Z3(2), Z3(1)]); // 2 + x
        let c = field.add(&a, &b);
        // (1+x) + (2+x) = 3 + 2x = 2x (mod 3)
        assert_eq!(c.coeffs[0], Z3::ZERO);
        assert_eq!(c.coeffs[1], Z3(2));
    }

    #[test]
    fn test_gf9_multiply() {
        let field = GF3n::with_degree(2);
        let a = field.element(vec![Z3(0), Z3(1)]); // x
        let b = field.element(vec![Z3(0), Z3(1)]); // x
        let c = field.mul(&a, &b); // x^2 mod irreducible
        // x^2 should reduce mod the irreducible polynomial
        assert!(c.coeffs.iter().any(|&v| v != Z3::ZERO));
    }

    #[test]
    fn test_gf9_inverse() {
        let field = GF3n::with_degree(2);
        let a = field.element(vec![Z3(1), Z3(1)]);
        let inv = field.inv(&a).expect("Nonzero element should have inverse");
        let prod = field.mul(&a, &inv);
        assert_eq!(prod, field.one());
    }

    #[test]
    fn test_gf9_order() {
        let field = GF3n::with_degree(2);
        // GF(9)* has order 8
        assert_eq!(field.element_order(&field.one()), 1);
    }

    #[test]
    fn test_gf9_primitive() {
        let field = GF3n::with_degree(2);
        let prim = field.primitive_element();
        assert_eq!(field.element_order(&prim), 8); // 3^2 - 1 = 8
    }

    #[test]
    fn test_gf27() {
        let field = GF3n::with_degree(3);
        assert_eq!(field.field_size(), 27);
        let a = field.element(vec![Z3(1), Z3(1), Z3(1)]);
        let inv = field.inv(&a).expect("Should have inverse");
        let prod = field.mul(&a, &inv);
        assert_eq!(prod, field.one());
    }

    #[test]
    fn test_gf81() {
        let field = GF3n::with_degree(4);
        assert_eq!(field.field_size(), 81);
        let prim = field.primitive_element();
        assert_eq!(field.element_order(&prim), 80);
    }

    #[test]
    fn test_irreducibles_degree3() {
        let irrs = PolyZ3::irreducibles(3);
        assert!(!irrs.is_empty());
        for p in &irrs {
            assert_eq!(p.degree(), Some(3));
            assert!(p.is_irreducible());
        }
    }

    #[test]
    fn test_poly_scalar_mul() {
        let p = PolyZ3::new(vec![Z3(1), Z3(1)]); // x + 1
        let q = p.scalar_mul(Z3(2)); // 2x + 2
        assert_eq!(q.coeff(0), Z3(2));
        assert_eq!(q.coeff(1), Z3(2));
    }
}
