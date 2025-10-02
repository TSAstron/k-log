use std::time::Instant;
//use rug::{Integer};//, Float, Complete};

use k_log::math_consts;
use k_log::{Frac, KFrac};
use k_log::KLog;

fn main() {
    let eta0 = Instant::now();
    let mut w = KFrac::from(math_consts::pi_ab);
    //let mut w = KFrac::from(math_consts::e_ab);
    //let mut w = KFrac::from(math_consts::zeta2_ab);
    for _ in 0..10_000 {
        w.prop();
    }
    w.gcd_div();
    w.stats();
    println!("{}", w.frac2.float(64));
    let mut x = KLog::from( w );
    while x.red1(0) {}
    let mut n0 = x.remainder.n;
    let mut g = n0 as f64 /x.tape.len() as f64;
    println!("tape len {} initial ratio {:.4}", x.tape.len(), g);
    let chunk = 10_000;
    for _ in 0..99 {
        x.red( chunk.min(10_000_000-x.tape.len()), Some(g), false);
        x.remainder.stats();
        g +=  (x.remainder.n-n0) as f64 /chunk as f64;
        g *= 0.5;
        println!("tape len {} convergent #{} ratio {:.4}", x.tape.len(), x.remainder.n, g);
        n0 = x.remainder.n;
        if x.tape.len() == 10_000_000 {break;}
    }
    let dl = x.remainder.frac2.den.to_string_radix(10).len();
    let t2 = eta0.elapsed().as_millis();
    println!("Klog has {} digits, convergent #{}, denominator len {}. {}ms", x.tape.len(), x.remainder.n, dl, t2 );
    println!("{}", x );
    let gt = x.gen_tape();
    println!("{:?}\n{} gigits, max: {:?}", &gt[99_900..100_010], gt.len(), gt.iter().max().unwrap() );
}
