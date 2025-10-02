// Implementation of Continued Logarithms for ordinary and continued fractions

use rug::Integer;

pub mod math_consts;
mod fractions;
mod series;
mod logs;

pub use fractions::{Frac, KFrac, RegFrac, RegFracOps};
pub use logs::{KLog, KLogOps};
pub use series::{Series};

impl KLogOps for Frac {
    fn iota( &mut self) {
        self.num -= &self.den;
        std::mem::swap( &mut self.den, &mut self.num);
    }
    fn theta( &mut self) {
        if self.num.is_even()  {
            self.num >>= 1;
            //self.num.div_exact_u_mut(2); // This is another nice option!
        } else {
            self.den <<= 1;
        }
    }
    fn gcd_div(&mut self) {
        self.gcd_div();
    }
}

impl KLog<Frac> {
    pub fn from_i128( a: i128, b: i128 ) -> Self {
        KLog::from( Frac::from_i128(a,b) )
    }
    pub fn from_int( a: Integer, b: Integer) -> Self {
        KLog::from( Frac::from_int(a,b) )
    }

    pub fn red1( &mut self) {
        if self.remainder.atleast2() {
            self.remainder.theta();
            self.tape.push(1);
        } else if self.remainder.atleast1() {
            self.remainder.iota();
            self.tape.push(0);
        } else {
            println!("\u{274c} Number fell below 1 during reduction!");
            panic!();
        }
    }
    pub fn red( &mut self, verbose: bool ) {
        if !self.remainder.atleast1() {
            println!("\u{274c} Can't reduce numbers (Frac) less than 1!");
            return;
        }
        while self.remainder.finite() {
            self.red1();
            if verbose { println!("{}", self); }
        }
        self.done = true;
    }

    fn reconstruct( mut x: Frac, s: &str ) -> Option<Frac> {
        for c in s.chars().rev() {
            if c == '0' {
                x.inv();
                x.num += &x.den;
            } else if c == '1' {
                x.num <<= 1;
            } else {
                println!("\u{274c} Unrecognized symbol '{}'", c);
                return None;
            }
        }
        x.gcd_div();
        Some(x)
    }

    pub fn parse( mut input: String) -> Option<(Frac,Frac)> {
    // A simple parser from k-log to ordinary fraction, gives Option
    // of the interval, possibly (x,x) if the number is fully determined
        match input.pop() {
            None => {
                println!("Empty k-log.");
            },
            Some('0') => {
                if let Some(rec) = Self::reconstruct( Frac::from(1), &input ) {
                    return Some((rec.clone(), rec));
                }
            },
            Some('1') => {
                println!("Last symbol '1'. Assuming incomplete expansion to produce an interval.");
                let x1 = Frac::from_i128(1,0);
                let x2 = Frac::from_i128(2,1);
                if let Some(y1) = Self::reconstruct( x1, &input ) {
                    if let Some(y2) = Self::reconstruct( x2, &input) {
                        return Some((y1, y2));
                    }
                }
            },
            Some(c) => {
                println!("\u{274c} Unrecognized symbol '{}' used.", c);
            }
        }
        None
    }
}

impl KLogOps for KFrac {
    fn iota( &mut self) {
        self.frac1.num -= &self.frac1.den;
        self.frac2.num -= &self.frac2.den;
        self.frac1.inv();
        self.frac2.inv();
    }
    fn theta( &mut self ) {
        if self.frac1.num.is_even() && self.frac2.num.is_even() {
            self.frac1.num >>= 1;
            self.frac2.num >>= 1;
        } else {
            self.frac1.den <<= 1;
            self.frac2.den <<= 1;
        }
    }
    fn gcd_div(&mut self) {
        self.gcd_div();
    }
    // Direct computation of the shorthand exapnsion eliminates all the powers of two in one go
    // but it requires the precomputation of s = log_2(x), which requires its own loop anyway...
    /*fn theta_iota( &mut self, mut s: u32 ) {
        let mut z = std::cmp::min( self.frac1.num.find_one(0).unwrap(), self.frac2.num.find_one(0).unwrap() );
        z = std::cmp::min( z, s );
        if z > 0 {
            self.frac1.num >>= z;
            self.frac2.num >>= z;
        }
        s -= z;
        if s > 0 {
            self.frac1.den <<= s;
            self.frac2.den <<= s;
        }
        self.iota();
    }*/
}

impl KLog<KFrac> {
    pub fn red1( &mut self, ratio: usize) -> bool {
        // This is the unchecked version, assuming propagation eventually helps
        // Infinite loop for numbers astronomically close to 1 -> [1-ε, 1+ε]
        let (s1, s2) = self.remainder.atleast2();
        if s1 && s2 {
            self.remainder.theta();
            self.tape.push(1);
            true
        } else if !s1 && !s2 && self.remainder.atleast1() == (true, true) {
            self.remainder.iota();
            self.tape.push(0);
            true
        } else {
            for _ in 0..ratio {
                self.remainder.prop();
            }
            false
        }
    }
    // irrationals require the target number of (new) digits to eventually stop
    pub fn red( &mut self, mut new_digits: usize, ratio: Option<f64>, verbose: bool ) {
        let limit = self.tape.len() + new_digits;
        if self.remainder.atleast1() == (false, false) {
            println!("\u{274c} KLog<KFrac>: Can't reduce numbers (KFrac) less than 1!");
            return;
        }
        if !self.remainder.finite() {
            println!("\u{274c} KLog<KFrac>: Both convergents are infinite!");
            return;
        }
        // pro-phase: extract what's freely available (but not more!)
        while new_digits > 0 && self.red1(0) {
            new_digits -= 1;
        }
        // pre-phase: assuming 1 term per 1 convergent by default 
        let rx = ratio.unwrap_or(1.0);
        for _ in 0..((rx * new_digits as f64).round() as usize) {
            self.remainder.prop();
        }
        self.remainder.gcd_div();
        let ry = rx.round().max(1.0) as usize ;
        if verbose {
            while self.tape.len() < limit {
                self.red1( ry );
                println!("{}", self );
            }
        } else {
            while self.tape.len() < limit {
                self.red1( ry );
            }
        }
    }
}
