fn main() {
    // Immutable variable (cannot change)
    let x = 5;
    println!("x = {}", x);

    // Mutable variable (can change)
    let mut counter = 0;
    println!("counter (start) = {}", counter);

    counter = counter + 1;
    println!("counter (after +1) = {}", counter);

    // With explicit type annotation
    let pi: f64 = 3.1415926535897932384626;
    println!("pi = {}", pi);
    let mut led_count: f64 = 6.000;
let mut voltage: f64 = 34.0000;
let mut is_on: bool = false;
println!("led_count = {}, voltage = {}, is_on = {}", led_count, voltage, is_on);
}
