// Generalized (or not) continued fractions for some mathematical constants

type I = i128;
type I2 = (i128, i128);

pub fn sqrt2_ab(n: I) -> I2 {
    if n == 0 {
        ( 1, 1 )
    } else {
        ( 1, 2 )
    }
}

pub fn sqrt13_ab(n: I) -> I2 {
    if n == 0 {
        (1, 3)
    } else if n%5 == 4 {
        (1, 6)
    } else {
        (1, 1)
    }
}

#[allow(unused_variables)]
pub fn gold_ab(n: I) -> I2 {
    (1, 1)
}

// The standard regular fraction
pub fn e_ab(n: I) -> I2 {
    if n == 0 {
        ( 1, 2 )
    } else if n%3 == 0 || n%3 == 1 {
        ( 1, 1 )
    } else {
        ( 1, 2*(n/3+1) )
    }
}

// A generalized continued fraction of e
pub fn e2_ab(n: I) -> I2 {
    if n == 0 {
        ( 1, 2 )
    } else if n == 1 {
        ( 1, 1 )
    } else {
        ( n-1, n )
    }
}

pub fn pi_ab(n: I) -> I2 {
    if n == 0 {
        ( 1, 0 )
    } else if n == 1 {
        (4, 1)
    } else {
        let m = n - 1;
        ( m*m, 2*m+1 )
    }
}

pub fn pi2_ab(n: I) -> I2 {
    if n == 0 {
        (1, 0)
    } else if n == 1 {
        (1008, -168)
    } else {
        (-1100736 + n*(4054680 + n*(-5999178 + n*(4556445 + n*(-1873686 + (396225 - 33750*n)*n)))),
         -3264 + n*(6008 + n*(-3637 + 725*n)) )
    }
}

// From ln(z) = 2 arth( (z-1)/(z+1) )
pub fn log2_ab(n: I) -> I2 {
    if n == 0 {
        (1, 0)
    } else if n == 1 {
        (2, 3)
    } else {
        (-(n-1)*(n-1), 6*n-3)
    }
}

// Accelerated by Apery, zeta(2) = pi^2/6
pub fn zeta2_ab(n: I) -> I2 {
    if n == 0 {
        (1, 0)
    } else if n == 1 {
        (5, 3)
    } else {
        let m = (n-1)*(n-1);
        ( m*m, 11*n*(n-1)+3 )
    }
}

// Accelerated by Apery, zeta(3)
pub fn zeta3_ab(n: I) -> I2 {
    if n == 0 {
        (1, 0)
    } else if n == 1 {
        (6, 5)
    } else {
        let m = (n-1)*(n-1);
        ( -m*m*m, (2*n-1)*(17*n*(n-1)+5) )
    }
}
