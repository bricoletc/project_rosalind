## Description 
Solutions to problems on [project Rosalind](https://rosalind.info/problems/list-view/)

## Structure

- Each subfolder is named according to the problem ID on the platform.
- Each subfolder contains an example dataset and a test (the one you download) dataset.
  There can be more than one of the latter if I failed preceding attempts!
- Each subfolder contains folders with code for solutions, named by language: e.g. `gasm/rust_gasm`

## To re-run

- For rust-based solutions:
  ```sh
  cd <problem_ID>/rust_<problem_ID>
  cargo test
  cargo run ../rosalind_<problem_ID>
  ```

