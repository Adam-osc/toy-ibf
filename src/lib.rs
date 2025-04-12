use fixedbitset::FixedBitSet;
use pyo3::exceptions::{PyKeyError, PyValueError};
use pyo3::prelude::*;
use rapidhash::RapidInlineHasher;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};

const GOLDEN_RATIO_64: u64 = 0x9e3779b97f4a7c15;

fn det_minimizer<'a>(
    window: &'a [u8],
    w: usize,
    k: usize,
    h_fn: &mut RapidInlineHasher,
    previous_minimizer_record: (&'a [u8], usize, usize),
) -> (&'a [u8], usize, usize) {
    let (mut prev_minimizer, mut prev_value, mut lifetime) = previous_minimizer_record;
    let mut lower_bound = 0;
    if lifetime > 0 {
        lifetime -= 1;
        lower_bound = w - 1;
    } else {
        prev_value = usize::MAX;
    }

    for i in lower_bound..w {
        let end_position = (i + k).min(window.len());
        let k_mer = &window[i..end_position];

        k_mer.hash(h_fn);
        let value = (h_fn.finish() ^ GOLDEN_RATIO_64) as usize;

        if value <= prev_value {
            prev_minimizer = k_mer;
            prev_value = value;
            lifetime = i;
        }

        if end_position == window.len() {
            break;
        }
    }

    (prev_minimizer, prev_value, lifetime)
}

fn hash(
    h_fn_1: &mut RapidInlineHasher,
    h_fn_2: &mut RapidInlineHasher,
    hashes: &mut Vec<usize>,
    fragment: &[u8],
) {
    fragment.hash(h_fn_1);
    let h1 = h_fn_1.finish();
    hashes.push(h1 as usize);

    fragment.hash(h_fn_2);
    let h2 = h_fn_2.finish();
    hashes.push(h2 as usize);

    for _ in 2..hashes.len() {
        hashes.push(h1.wrapping_add(h2).rotate_left(5) as usize);
    }
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|c| if c & 2 != 0 { c ^ 4 } else { c ^ 21 })
        .collect()
}

#[pyclass]
struct InterleavedBloomFilter {
    fragment_length: usize,
    overlap: usize,
    w: usize,
    k: usize,
    num_of_bins: usize,
    single_filter_size: usize,
    preserved_pct: f64,
    max_idx: usize,
    bin_idxs: HashMap<String, usize>,
    active_filter: Vec<FixedBitSet>,
    filters: HashMap<(String, String), Vec<FixedBitSet>>,
    h_fn_1: RapidInlineHasher,
    h_fn_2: RapidInlineHasher,
    hashes: Vec<usize>,
    counting_vector: Vec<usize>,
    binning_vector: FixedBitSet,
}

#[pymethods]
impl InterleavedBloomFilter {
    #[new]
    fn new(
        num_of_bins: usize,
        single_filter_size: usize,
        fragment_length: usize,
        overlap: usize,
        w: usize,
        k: usize,
        num_hashes: usize,
        preserved_pct: f64
    ) -> PyResult<Self> {
        if num_hashes < 2 {
            return Err(PyValueError::new_err(
                "Atleast two hash functions have to be used",
            ));
        }
        if overlap >= fragment_length {
            return Err(PyValueError::new_err(
                "Overlap cannot be greater or equal to the fragment length",
            ));
        }
        if preserved_pct <= 0.0 || preserved_pct >= 1.0 {
            return Err(PyValueError::new_err(
                "Preserved percentage must be a probability between 0 and 1 (exclusive)",
            ));
        }
        // There is are no single bloom filters...
        // But the overall interleaved bloom filter size is: <single_filter_size> * <num_of_bins>
        let max_idx = 0;
        let active_filter: Vec<FixedBitSet> = (0..single_filter_size)
            .map(|_| {
                let mut bit_set = FixedBitSet::with_capacity(num_of_bins);
                bit_set.set_range(.., false);
                bit_set
            })
            .collect();
        let filters: HashMap<(String, String), Vec<FixedBitSet>> = HashMap::new();
        let bin_idxs: HashMap<String, usize> = HashMap::new();

        let hashes: Vec<usize> = Vec::with_capacity(num_hashes);
        let h_fn_1 = RapidInlineHasher::default();
        let h_fn_2 = RapidInlineHasher::new(GOLDEN_RATIO_64);

        let counting_vector: Vec<usize> = vec![0; num_of_bins];
        let binning_vector = FixedBitSet::with_capacity(num_of_bins);

        Ok(Self {
            fragment_length,
            overlap,
            w,
            k,
            num_of_bins,
            single_filter_size,
            max_idx,
            bin_idxs,
            active_filter,
            preserved_pct,
            filters,
            h_fn_1,
            h_fn_2,
            hashes,
            counting_vector,
            binning_vector,
        })
    }
    pub fn insert_sequence(&mut self, seq_id: (String, String), seq: &str) {
        if self.filters.contains_key(&seq_id) {
            return;
        }
        let bin_idx = self.bin_idxs.entry(seq_id.0.clone()).or_insert(self.max_idx);
        *bin_idx += 1;
        
        let mut filter: Vec<FixedBitSet> = (0..self.single_filter_size)
            .map(|_| {
                let mut bit_set = FixedBitSet::with_capacity(self.num_of_bins);
                bit_set.set_range(.., false);
                bit_set
            })
            .collect();

        for frag_start in (0..seq.len())
            .step_by(self.fragment_length - self.overlap)
        {
            let fragment = &seq[frag_start..(frag_start + self.fragment_length).min(seq.len())];
            let num_windows = if fragment.len() >= (self.w + self.k - 1) {
                fragment.len() - (self.w + self.k - 1)
            } else {
                0
            } + 1;

            let mut prev_minimizer: &[u8] = &[];
            let mut prev_value: usize = usize::MAX;
            let mut lifetime: usize = 0;

            for j in 0..num_windows {
                let window = fragment[j..(j + self.w + self.k - 1).min(fragment.len())].as_bytes();
                let (minimizer, value, new_lifetime) = det_minimizer(
                    window,
                    self.w,
                    self.k,
                    &mut self.h_fn_1,
                    (prev_minimizer, prev_value, lifetime),
                );

                if lifetime > 0 && prev_value == value && prev_minimizer == minimizer {
                    lifetime = new_lifetime;
                    continue;
                }
                prev_minimizer = minimizer;
                prev_value = value;
                lifetime = new_lifetime;

                for mer in [prev_minimizer, &reverse_complement(prev_minimizer)] {
                    hash(&mut self.h_fn_1, &mut self.h_fn_2, &mut self.hashes, mer);
                    for digest in self.hashes.drain(..) {
                        FixedBitSet::set(&mut filter[digest % self.single_filter_size], *bin_idx % self.num_of_bins, true);
                    }
                }
            }
        }

        self.filters.insert(seq_id, filter);
    }

    pub fn activate_filter(&mut self, seq_id: (String, String)) -> PyResult<()> {
        let filters = self.filters.get(&seq_id).ok_or_else(|| {
            PyKeyError::new_err("Tried to activate a filter that is not present.")
        })?;

        for (i, filter) in filters.iter().enumerate() {
            self.active_filter[i] |= filter;
        }
        Ok(())
    }

    pub fn is_sequence_present(&mut self, seq: &str) -> bool {
        let mut max_count: usize = 0;
        let num_windows = if seq.len() >= (self.w + self.k - 1) {
            seq.len() - (self.w + self.k - 1)
        } else {
            0
        } + 1;

        let threshold = (num_windows as f64 * self.preserved_pct).ceil() as usize;
        self.counting_vector.fill(0);

        let mut prev_minimizer: &[u8] = &[];
        let mut prev_value: usize = usize::MAX;
        let mut lifetime: usize = 0;

        for i in 0..num_windows {
            let window = &seq[i..(i + self.w + self.k - 1).min(seq.len())];
            let (minimizer, value, new_lifetime) = det_minimizer(
                window.as_bytes(),
                self.w,
                self.k,
                &mut self.h_fn_1,
                (prev_minimizer, prev_value, lifetime),
            );

            if lifetime > 0 && prev_value == value && prev_minimizer == minimizer {
                lifetime = new_lifetime;
            } else {
                prev_minimizer = minimizer;
                prev_value = value;
                lifetime = new_lifetime;

                self.binning_vector.set_range(.., true);
                hash(
                    &mut self.h_fn_1,
                    &mut self.h_fn_2,
                    &mut self.hashes,
                    prev_minimizer,
                );
                for digest in self.hashes.drain(..) {
                    self.binning_vector &= &self.active_filter[digest % self.single_filter_size];
                }
            }

            for j in self.binning_vector.ones() {
                self.counting_vector[j] += 1;
                max_count = max_count.max(self.counting_vector[j]);
                if max_count >= threshold {
                    return true;
                }
                if threshold - max_count > num_windows - i - 1 {
                    return false;
                }
            }
        }
        false
    }
}

#[pymodule]
fn interleaved_bloom_filter(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<InterleavedBloomFilter>()?;
    Ok(())
}
