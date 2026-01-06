// This binary is deprecated; the Python CLI (`scnado`) is the supported entry point.
// A stub main is kept so cargo-clippy and other tooling have a valid binary target.

fn main() {
    eprintln!("The Rust CLI is deprecated. Use the Python CLI: scnado <command> ...");
}
