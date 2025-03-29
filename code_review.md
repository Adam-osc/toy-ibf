# Code review
I went into details since the code is for your Bachelor thesis.
Also the suggested code was not tested, so look if the logic still works.

## Error handling
Currently uses `panic!` which abruptly terminates the program therefore being unsafe.
Could use `.expect()` or `Result` type for more graceful error management.

```rust
pub fn activate_filter(&mut self, seq_id: &str) {
    let filter = self.filters.get(seq_id)
        .expect(&format!("Tried to activate a filter that is not present. Filter '{}' not found.", seq_id));

    for (i, f) in filter.iter().enumerate() {
        self.active_filter[i] |= f;
    }
}
```

or

```rust
pub fn activate_filter(&mut self, seq_id: &str) -> Result<(), String> {
    self.filters.get(seq_id)
        .map(|filter| {
            for (i, f) in filter.iter().enumerate() {
                self.active_filter[i] |= f;
            }
            Ok(())
        })
        .unwrap_or_else(|| Err(format!("Tried to activate a filter that is not present. Filter '{}' not found.", seq_id)))
}
```

## Code duplication
The filter initialization happens in multiple parts of the code.
Could extract it into a separate method.

```rust
fn create_empty_filter(num_vectors: usize, size_per_vector: usize) -> Vec<FixedBitSet> {
    (0..num_vectors).map(|_| FixedBitSet::with_capacity(size_per_vector)).collect()
}

let mut filter = create_empty_filter(self.num_vectors, self.size_per_vector);
```

Also in `is_sequence_present`, `seq.as_bytes()` is repeatedly called.
Could save the result in a variable similarly to `threshould`.

```rust
pub fn is_sequence_present(&mut self, seq: &str) -> bool {
    let (_, upper) = calculate_confidence_intervals(self.standard_norm_dist, seq.len(), self.k, self.error_rate, self.confidence);
    let threshold = seq.len() - self.k + 1 - upper;
    let seq_bytes = seq.as_bytes();
    self.is_sequence_present_aux(seq.as_bytes(), threshold) || self.is_sequence_present_aux(&reverse_complement(seq.as_bytes()), threshold)
}
```

## Memory allocation
The hashes vector is drained repeatedly, which reallocates memory on each call.
Instead of `.drain(..)`, reuse the vector.

```rust
for h in &self.hashes {
    filter[h % self.num_vectors].set(i, true);
}

// Reset without reallocation
self.hashes.clear();
```

## Documentation and naming
Lacks comprehensive documentation comments.
No clear explanations of method purposes and parameter meanings.
Could benefit from comments explaining algorithm specifics.

Some variable names make it hard to understand their use (`h_fn_1`, `h_fn_2`).
Could use more descriptive names like `hash_fn_1`, `hash_fn_2` to add a bit more clarity.

### Note: Commenting in Rust

#### 1. `//` - Regular Comments
- Used for inline comments that do not generate documentation.
- Meant for developers reading the source code.

#### 2. `///` - Documentation Comments (Outer Doc)
- Used to document functions, structs, modules, etc.
- Supports Markdown formatting.
- Included in generated documentation when using `cargo doc`.

## Potential optimizations suggested by AI
- There is a lot of cases where the type is not defined which can be unsafe. Explicit typing is safer.
- Consider using SmallVec for small hash collections to avoid heap allocations
- Explore using SIMD instructions for faster bit operations
- Implement parallel processing for large sequence sets
