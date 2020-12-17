
// Some helper math
// C-style call syntax, e.g. `min(a,b)` instead of `a.min(b)`
pub(crate) fn sqrt(a: f64) -> f64 {
    a.sqrt()
}
pub(crate) fn log(a: f64) -> f64 {
    a.ln()
}
pub(crate) fn exp(a: f64) -> f64 {
    a.exp()
}
pub(crate) fn MAX(a: f64, b: f64) -> f64 {
    a.max(b)
}
pub(crate) fn MIN(a: f64, b: f64) -> f64 {
    a.min(b)
}
pub(crate) fn abs(a: f64) -> f64 {
    a.abs()
}
pub(crate) fn atan(a: f64) -> f64 {
    a.atan()
}
pub(crate) fn pow(a: f64, b: f64) -> f64 {
    a.powf(b)
}
pub(crate) fn max(a: f64, b: f64) -> f64 {
    a.max(b)
}
