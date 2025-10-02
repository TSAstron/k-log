// Ordinary and continued fractions

use rug::{Integer, Float, Complete};
use std::fmt::Display;

#[derive(Clone, Debug, PartialEq)]
pub struct Frac {
    pub num: Integer,
    pub den: Integer,
}

impl Frac {
    pub fn from_i128(x: i128, y: i128) -> Self {
        Self {num: Integer::from(x), den: Integer::from(y)}
    }
    pub fn from_int(x: Integer, y: Integer) -> Self {
        Self {num: x, den: y}
    }
    pub fn float (&self, prec: u32) -> Float {
        Float::with_val(prec, &self.num)/&self.den 
    }
    pub fn inv( &mut self ) {
        std::mem::swap( &mut self.num, &mut self.den );
    }
    pub fn gcd_div( &mut self ) {
        // reduce the fraction via gcd
        let g = self.num.gcd_ref( &self.den ).complete();
        self.num /= &g;
        self.den /= &g;
    }
    pub fn atleast1( &self ) -> bool {
        if self.den >= 0 { // We still want inf >= 1
            self.den >= 0 && self.num >= self.den
        } else {
            self.den < 0 && self.num <= self.den
        }
    }
    pub fn atleast2( &self ) -> bool {
        if self.den >= 0 { // We still want inf >= 2
            Integer::from(&self.num - &self.den) >= self.den
        } else {
            Integer::from(&self.num - &self.den) <= self.den
        }
    }
    pub fn finite( &self ) -> bool {
        self.den != 0
    }
}

impl From<i128> for Frac {
    fn from( a: i128 ) -> Self {
        Frac::from_i128(a, 1)
    }
}

impl From<Integer> for Frac {
    fn from( a: Integer ) -> Self {
        Frac::from_int(a, Integer::from(1))
    }
}

impl Display for Frac {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let mut num_str = self.num.to_string();
        let mut den_str = self.den.to_string();
        if num_str.len() > 64 {
            num_str = format!("{}...{}", &num_str[..32], &num_str[num_str.len()-32..]);
        }
        if den_str.len() > 64 {
            den_str = format!("{}...{}", &den_str[..32], &den_str[den_str.len()-32..]);
        }
        write!( f, "{}/{}", num_str, den_str )
    }
}

type PropF = fn(i128) -> (i128, i128);
type Mat<T> = ((T, T), (T, T));

pub struct KFrac {
    pub n: i128,
    ab: PropF,
    pub frac1: Frac,
    pub frac2: Frac,
}

impl Display for KFrac {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        // Proper order but only for actual convergents, KLog reductions spoil it
        if self.n % 2 == 0 {
            write!( f, "[{}; {}]", self.frac2, self.frac1 )
        } else {
            write!( f, "[{}; {}]", self.frac1, self.frac2 )
        }
    }
}

impl KFrac {
    pub fn from( ab: PropF ) -> Self {
        KFrac {
            n: 0,
            ab: ab,
            // New fractions are always: | 1  b_0 |
            //                           | 0  1   |
            frac1: Frac::from_i128( 1, 0),
            frac2: Frac::from_i128( ab(0).1 , 1),
        }
    }
    pub fn atleast1(&self) -> (bool, bool) {
        (self.frac1.atleast1() , self.frac2.atleast1() )
    }
    pub fn atleast2(&self) -> (bool, bool) {
        (self.frac1.atleast2() , self.frac2.atleast2() )
    }
    pub fn finite(&self) -> bool {
    // Only when both are 1/0 can we be sure it's infinite
        self.frac1.finite() || self.frac2.finite()
    }
    pub fn gcd_div(&mut self) {
        // divide both convergents by gcd of all 4 components
        let g1 = self.frac1.num.gcd_ref( &self.frac1.den ).complete();
        let mut g = self.frac2.num.gcd_ref( &self.frac2.den ).complete();
        g.gcd_mut( &g1 );
        self.frac1.num /= &g;
        self.frac1.den /= &g;
        self.frac2.num /= &g;
        self.frac2.den /= &g;
    }

    pub fn prop_bs(&mut self, new: i128) {
        // binary splitting propagation, to advance 'self' by 'new' convergents
        fn mat_mul( ((a11, a12), (a21, a22)): Mat<Integer>,
                    ((b11, b12), (b21, b22)): Mat<Integer> ) -> Mat<Integer> {
            let s11 = a11.clone() * &b11 + a12.clone() * &b21;
            let s12 = (a11 * &b12) + (a12 * &b22);
            let s21 = a21.clone() * b11 + a22.clone() * b21;
            let s22 = a21 * b12 + a22 * b22;
           ( (s11, s12), (s21, s22) ) 
        }

        fn matrix_bs(f: PropF, m1: i128, m2: i128) -> Mat<Integer> {
            if m1 == m2 {
                let (a,b) = f(m1);
                ( (Integer::from(0), Integer::from(a)), 
                  (Integer::from(1), Integer::from(b)) )
            } else {
                let mid = m1 + (m2-m1)/2;
                mat_mul(matrix_bs(f, m1, mid), matrix_bs(f, mid + 1, m2))
            }
        }

        ((self.frac1.num, self.frac2.num), 
         (self.frac1.den, self.frac2.den)) = mat_mul( 
            ((self.frac1.num.clone(), self.frac2.num.clone()), 
             (self.frac1.den.clone(), self.frac2.den.clone())), 
            matrix_bs(self.ab, self.n + 1, self.n + new) );

        //self.frac1.num = a11;
        //self.frac2.num = a12;
        //self.frac1.den = a21;
        //self.frac2.den = a22;
        self.n += new;
    }

    pub fn prop(&mut self) {
        self.n += 1;
        let (a, b) = (self.ab)(self.n);
        self.frac1.num *= &a;
        self.frac1.den *= &a;
        self.frac1.num += &b * &self.frac2.num;
        self.frac1.den += &b * &self.frac2.den;
        std::mem::swap( &mut self.frac1, &mut self.frac2 );
    }

    pub fn stats(&self) -> (u32, f64) {
        // Some rudimentary error and precision estimates
        let den_len = &self.frac2.den.to_string_radix(10).len();
        let err_den = Integer::from( &self.frac2.den * &self.frac1.den );
        if err_den == 0 {
            println!( "[convergent {}] err around ? (inf), denominator of length {}", self.n, den_len );
            return (0, 0.0)
        }
        let err_num = (Integer::from(&self.frac2.num * &self.frac1.den) - Integer::from(&self.frac1.num * &self.frac2.den)).abs();
        let whole = Integer::from( &self.frac2.num / &self.frac2.den ).significant_bits();
        let mut prec2 = whole + err_den.significant_bits() + 1;
        prec2 = match prec2.checked_sub( err_num.significant_bits() ) {
            Some(z) => z,
            None => 0 // the rare case when the error is larger than the number itself (frac2)
        };
        let err = Float::with_val(prec2+3, err_num)/err_den;
        let decims = -err.clone().log10().to_f64().floor();
        println!( "[convergent {}] err ~ decimal place {} ({:.2e}), denominator of length {}", self.n, decims, err, den_len );
        ( prec2, decims )
    }
}

pub trait RegFracOps {}

pub struct RegFrac<T: RegFracOps> {
    pub tape: Vec<u128>,
    pub remainder: T,
}

impl<T: RegFracOps> From<T> for RegFrac<T> {
    fn from( a: T ) -> Self {
        RegFrac { tape: Vec::new(), remainder: a }
    }
}

impl<T> Display for RegFrac<T>
where T: RegFracOps + Display {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!( f, "{:?} ({})", self.tape, self.remainder)
    }
}

impl RegFracOps for Frac {}

impl RegFrac<Frac> {
    pub fn red1(&mut self) {
        let rem = &mut self.remainder;
        let (a, b) = rem.num.div_rem_ref(&rem.den).complete();
        if let Some(c) = a.to_u128() {
            self.tape.push(c);
            rem.num = b;
            std::mem::swap( &mut rem.num, &mut rem.den);
        } else {
            println!("The next term does not fit in u128.");
            panic!();
        }
    }

    pub fn red(&mut self, verbose: bool) {
        if verbose {
            while self.remainder.finite() {
                self.red1();
                println!("{}", self);
            }
        } else {
            while self.remainder.finite() {
                self.red1();
            }
        }
    }
}

impl RegFracOps for KFrac {}

impl RegFrac<KFrac> {
    pub fn red1(&mut self, ratio: usize) -> bool {
        let rem = &mut self.remainder;
        if !rem.frac1.finite() || !rem.frac2.finite() {
            rem.prop();
            return false;
        }
        let (a1, b1) = rem.frac1.num.div_rem_ref(&rem.frac1.den).complete();
        let (a2, b2) = rem.frac2.num.div_rem_ref(&rem.frac2.den).complete();
        if a1 == a2 {
            if let Some(c) = a1.to_u128() {
                self.tape.push(c);
                rem.frac1.num = b1;
                std::mem::swap( &mut rem.frac1.num, &mut rem.frac1.den);
                rem.frac2.num = b2;
                std::mem::swap( &mut rem.frac2.num, &mut rem.frac2.den);
                return true;
            } else {
                println!("The next term does not fit in u128.");
                panic!();
            }
        } else {
            for _ in 0..ratio {
                rem.prop();
            }
            return false;
        }
    }

    pub fn red_adapt(&mut self, new_terms: usize, chunk: usize) {
        if new_terms <= chunk { // 10_000 seems like a good choice
            self.red( new_terms, None )
        } else {
            let mut l0 = self.tape.len();
            let limit = l0 + new_terms;
            let mut n0 = self.remainder.n;
            for _ in 0..chunk {
                self.remainder.prop()
            }
            self.remainder.gcd_div(); // Move it Afterwards???
            while self.red1(0) {} // Extract as many as possible without prop
            let mut l1 = self.tape.len();
            let mut g = (self.remainder.n - n0) as f32 / (l1 - l0 ) as f32;
            println!("{}", g);
            while limit > l1 {
                n0 = self.remainder.n;
                l0 = l1;
                let loc_limit = (limit - l1).min(chunk) as f32 * g;
                for _ in 0..(loc_limit as usize) {
                    self.remainder.prop();
                }
                self.remainder.gcd_div();
                let g0 = (g as usize).max(1);
                for _ in 0..(limit-l1) {
                    if !self.red1( g0 ) {
                        break;
                    }
                }
                l1 = self.tape.len();
                g = (self.remainder.n - n0) as f32 / (l1-l0) as f32;
                println!("tape len {} converget #{} ratio {}", l1, self.remainder.n, g);
            }
        }
    }
    pub fn red(&mut self, new_terms: usize, ratio: Option<f32> ) {
        let limit = self.tape.len() + new_terms;
        // assuming 1 KFrac step produces ~1 regular fraction term is a reasonable default
        let r0 = ratio.unwrap_or(1.0);
        let scaled = (r0 * new_terms as f32).ceil() as usize;
        for _ in 0..scaled {
            self.remainder.prop();
        }
        self.remainder.gcd_div();
        // prop is repeated rx times inside red1
        let rx  = r0.ceil().max(1.0) as usize;
        while self.tape.len() < limit {
            self.red1(rx);
        }
    }
}
