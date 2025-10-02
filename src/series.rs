// Series with rational terms

use crate::fractions::Frac;

type PropS = fn(i128, &Frac) -> Frac;
type PropSm = fn(i128, &mut Frac);

#[derive(Debug)]
pub struct Series {
    pub n: i128,
    u: PropSm,
    pub curr: Frac,
    pub err: PropS,
    pub sum: Frac,
}

impl Series {
    pub fn from( u: PropSm, err: PropS ) -> Self {
        let mut undef = Frac::from(0);
        u(0, &mut undef); // u_0 should not depend on the value of undef (u_{-1})
        Series { 
            n: 0, 
            u: u, 
            err: err, 
            sum: undef.clone(),
            curr: undef }
    }

    pub fn error( &self ) -> Frac {
        (self.err)( self.n, &self.curr )
    }

    pub fn prop( &mut self ) {
        self.n += 1;
        (self.u)(self.n, &mut self.curr);
        self.sum.num *= &self.curr.den;
        self.sum.num += &self.curr.num * &self.sum.den;
        self.sum.den *= &self.curr.den;
    }
}
