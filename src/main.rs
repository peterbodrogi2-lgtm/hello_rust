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
}