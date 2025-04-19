use fixedbitset::FixedBitSet;
use pyo3::exceptions::{PyKeyError,PyIndexError, PyValueError};
use pyo3::prelude::*;
use rapidhash::{rapidhash, rapidhash_seeded};
use std::collections::HashMap;

const GOLDEN_RATIO_64: u64 = 0x9e3779b97f4a7c15;

fn det_minimizer<'a>(
    window: &'a [u8],
    rc_window: &'a [u8],
    w: usize,
    k: usize,
    previous_minimizer_record: (&'a [u8], usize, usize)
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

        for k_mer in [&window[i..end_position], &rc_window[i..end_position]] {
            let value = (rapidhash(k_mer) ^ GOLDEN_RATIO_64) as usize;

            if value <= prev_value {
                prev_minimizer = k_mer;
                prev_value = value;
                lifetime = i;
            }
        }

        if end_position == window.len() {
            break;
        }
    }

    (prev_minimizer, prev_value, lifetime)
}

fn hash(
    hashes: &mut Vec<usize>,
    fragment: &[u8],
) {
    let h1 = rapidhash(fragment);
    hashes.push(h1 as usize);
    let h2 = rapidhash_seeded(fragment, GOLDEN_RATIO_64);
    hashes.push(h1 as usize);

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
    w: usize,
    k: usize,
    num_of_bins: usize,
    single_filter_size: usize,
    max_idx: usize,
    bin_idxs: HashMap<String, usize>,
    bin_names: Vec<String>,
    active_filter: Vec<FixedBitSet>,
    filters: HashMap<String, Vec<FixedBitSet>>,
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
        w: usize,
        k: usize,
        num_hashes: usize
    ) -> PyResult<Self> {
        if num_hashes < 2 {
            return Err(PyValueError::new_err(
                "Atleast two hash functions have to be used",
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
        let filters: HashMap<String, Vec<FixedBitSet>> = HashMap::new();
        let bin_idxs: HashMap<String, usize> = HashMap::new();
        let bin_names: Vec<String> = Vec::with_capacity(num_of_bins);

        let hashes: Vec<usize> = Vec::with_capacity(num_hashes);

        let counting_vector: Vec<usize> = vec![0; num_of_bins];
        let binning_vector = FixedBitSet::with_capacity(num_of_bins);

        Ok(Self {
            fragment_length,
            w,
            k,
            num_of_bins,
            single_filter_size,
            max_idx,
            bin_idxs,
            bin_names,
            active_filter,
            filters,
            hashes,
            counting_vector,
            binning_vector,
        })
    }
    pub fn insert_sequence(&mut self, bin_id: String, seq: &str) -> PyResult<()> {
        if self.max_idx >= self.num_of_bins {
            return Err(PyIndexError::new_err("Attempted to insert into a non‑existent bin"));
        }

        self.bin_names.push(bin_id.clone());
        let bin_idx = self.bin_idxs.entry(bin_id.clone()).or_insert_with(|| {
            let i = self.max_idx;
            self.max_idx += 1;
            i
        });

        let rc_seq = reverse_complement(seq.as_bytes());
        let mut filter: Vec<FixedBitSet> = (0..self.single_filter_size)
            .map(|_| {
                let mut bit_set = FixedBitSet::with_capacity(self.num_of_bins);
                bit_set.set_range(.., false);
                bit_set
            })
            .collect();

        for frag_start in (0..seq.len())
            .step_by(self.fragment_length)
        {
            let fragment = &seq[frag_start..(frag_start + self.fragment_length).min(seq.len())];
            let rc_fragment = &rc_seq[frag_start..(frag_start + self.fragment_length).min(rc_seq.len())];
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
                let rc_window = &rc_fragment[j..(j + self.w + self.k - 1).min(rc_fragment.len())];
                let (minimizer, value, new_lifetime) = det_minimizer(
                    window,
                    rc_window,
                    self.w,
                    self.k,
                    (prev_minimizer, prev_value, lifetime),
                );

                if lifetime > 0 && prev_value == value && prev_minimizer == minimizer {
                    lifetime = new_lifetime;
                    continue;
                }
                prev_minimizer = minimizer;
                prev_value = value;
                lifetime = new_lifetime;

                hash(&mut self.hashes, prev_minimizer);
                for digest in self.hashes.drain(..) {
                    FixedBitSet::set(&mut filter[digest % self.single_filter_size], *bin_idx, true);
                }
            }
        }

        self.filters.insert(bin_id, filter);
        Ok(())
    }

    pub fn activate_filter(&mut self, bin_id: String) -> PyResult<()> {
        let filters = self.filters.get(&bin_id).ok_or_else(|| {
            PyKeyError::new_err("Tried to activate a filter that is not present")
        })?;

        for (i, filter) in filters.iter().enumerate() {
            self.active_filter[i] |= filter;
        }
        Ok(())
    }

    pub fn reset_filter(&mut self) {
        for single_filter in self.active_filter.iter_mut() {
            single_filter.set_range(.., false);
        }
    }

    pub fn is_sequence_present(&mut self, seq: &str, preserved_pct: f64) -> PyResult<Option<&str>> {
        if preserved_pct <= 0.0 || preserved_pct >= 1.0 {
            return Err(PyValueError::new_err(
                "Preserved percentage must be a probability between 0 and 1 (exclusive)",
            ));
        }

        let mut max_count: usize = 0;
        let mut max_present_count: usize = 0;
        let mut max_present_bin_id: Option<&str> = None;
        let num_windows = if seq.len() >= (self.w + self.k - 1) {
            seq.len() - (self.w + self.k - 1)
        } else {
            0
        } + 1;

        let rc_seq = reverse_complement(seq.as_bytes());
        let threshold = (num_windows as f64 * preserved_pct).ceil() as usize;
        self.counting_vector.fill(0);

        let mut prev_minimizer: &[u8] = &[];
        let mut prev_value: usize = usize::MAX;
        let mut lifetime: usize = 0;

        for i in 0..num_windows {
            let window = &seq[i..(i + self.w + self.k - 1).min(seq.len())];
            let rc_window = &rc_seq[i..(i + self.w + self.k - 1).min(seq.len())];
            let (minimizer, value, new_lifetime) = det_minimizer(
                window.as_bytes(),
                rc_window,
                self.w,
                self.k,
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

                if max_count >= threshold && self.counting_vector[j] > max_present_count {
                    max_present_count = self.counting_vector[j];
                    max_present_bin_id = Some(&self.bin_names[j]);
                }

                if max_count < threshold && threshold - max_count > num_windows - i - 1 {
                    return Ok(None);
                }
            }
        }
        Ok(max_present_bin_id)
    }
}

#[pymodule]
fn interleaved_bloom_filter(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<InterleavedBloomFilter>()?;
    Ok(())
}
