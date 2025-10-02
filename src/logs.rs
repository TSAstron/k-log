// Continued logarithms - the struct, and trait specifying the required operations for a type to be
// usable as the remainder of a KLog

use std::fmt::Display;

pub trait KLogOps {
    // The operation corresponding to '0': x -> 1/(x-1)
    fn iota(&mut self);
    // The operation corresponding to '1': x -> x/2
    fn theta(&mut self);
    // reducing the remaining fraction(s) via gcd -- costly! Currently only used at creation.
    fn gcd_div(&mut self);
    //fn theta_iota(&mut self, &mut usize); // composition of the above
}

pub struct KLog<T: KLogOps> {
    pub tape: Vec<u8>,
    pub remainder: T,
    pub done: bool,
}

impl<T> Display for KLog<T>
where T: KLogOps + Display {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.tape.len() == 0 {
            return write!( f, "_ ({})", self.remainder );
        }
        let tape_s = self.tape_str();
        let tlen = tape_s.len();
        if tlen > 64 {
            write!( f, "{}...{} ({})", &tape_s[..32], &tape_s[tlen-32..], self.remainder )
        } else {
            write!( f, "{} ({})", tape_s, self.remainder )
        }
    }
}

impl<T: KLogOps> KLog<T> {
    pub fn tape_str( &self ) -> String {
    // k-log digits as String
        self.tape.iter().map( |q| q.to_string() ).collect()
    }

    pub fn pretty( &self ) -> String {
    // pretty expansion including the end symbol 
        let mut prettape = String::new();
        for c in self.tape.iter() {
            match c {
                0 => prettape.push('○'), // ❍
                1 => prettape.push('△'),
                _ => {println!("\u{274c} This shouldn't happen!");}
            }
        }
        if self.done {
            prettape.push('□'); // ∞
        }
        prettape
    }

    pub fn gen_tape( &self ) -> Vec<u128> {
    // k-log in shorthand notation
        let mut gen_tape = Vec::<u128>::new();
        let mut g = 0;
        for &s in self.tape.iter() {
            if s == 1 {
                g += 1;
            } else {
                gen_tape.push(g);
                g = 0;
            }
        }
        if g > 0 {
            gen_tape.push(g);
        }
        gen_tape
    }
}

impl<T: KLogOps> From<T> for KLog<T> {
    fn from( mut f: T ) -> Self {
        f.gcd_div();
        KLog { tape: Vec::new(), remainder: f, done: false }
    }
}
