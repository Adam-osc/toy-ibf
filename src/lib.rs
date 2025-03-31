use fixedbitset::FixedBitSet;
use pyo3::prelude::*;
use rapidhash::RapidInlineHasher;
use statrs::distribution::{ContinuousCDF, Normal};
use std::collections::HashMap;
use std::hash::{Hash, Hasher};

fn calc_num_fragments(seq_length: usize, fragment_length: usize, overlap: usize) -> usize {
    (seq_length / (fragment_length - overlap)) + (seq_length % (fragment_length - overlap) == 0) as usize
}

fn calc_size_per_vector(num_hashes: usize, fragment_length: usize, k: usize) -> usize {
    let r = 0.01_f64.powf(1.0 / num_hashes as f64);
    let max_k_mer_num = fragment_length - k + 1;
    (-1.0 / ((1.0 - r).powf(1.0 / (num_hashes * max_k_mer_num) as f64) - 1.0)).ceil() as usize
}

fn hash(h_fn_1: &mut RapidInlineHasher, h_fn_2: &mut RapidInlineHasher, hashes: &mut Vec<usize>, fragment: &[u8]) {
    fragment.hash(h_fn_1);
    let h1 =  h_fn_1.finish();
    hashes.push(h1 as usize);

    fragment.hash(h_fn_2);
    let h2 = h_fn_2.finish();
    hashes.push(h2 as usize);

    for _ in 2..hashes.len() {
        hashes.push(h1.wrapping_add(h2).rotate_left(5) as usize);
    }
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|c| if c & 2 != 0 { c ^ 4 } else { c ^ 21 })
        .collect()
}

fn calculate_confidence_intervals(standard_norm_dist: Normal, seq_length: usize, k: usize, error_rate: f64, confidence: f64) -> (usize, usize) {
    let q = 1.0 - (1.0 - error_rate).powi(k as i32);
    let num_k_mers = seq_length - k + 1;
    let variance = num_k_mers as f64 * (1.0 - q) * (q * (2.0 * k as f64 + (2.0 / error_rate) - 1.0) - 2.0 * k as f64)
        + (k * (k - 1)) as f64 * (1.0 - q).powi(2)
        + ((2.0 * (1.0 - q)) / error_rate.powi(2)) * ((1.0 + (k - 1) as f64 * (1.0 - q)) * error_rate - q);
    let alpha = 1.0 - confidence;
    let z = standard_norm_dist.inverse_cdf(1.0 - alpha / 2.0);
    let mean = num_k_mers as f64 * q;
    ((mean - z * variance.sqrt()).floor() as usize, (mean + z * variance.sqrt()).ceil() as usize)
}

#[pyclass]
struct InterleavedBloomFilter {
    fragment_length: usize,
    overlap: usize,
    k: usize,
    num_vectors: usize,
    size_per_vector: usize,
    error_rate: f64,
    confidence: f64,
    active_filter: Vec<FixedBitSet>,
    filters: HashMap<(String, String), Vec<FixedBitSet>>,
    h_fn_1: RapidInlineHasher,
    h_fn_2: RapidInlineHasher,
    hashes: Vec<usize>,
    standard_norm_dist: Normal,
    counting_vector: Vec<usize>,
    binning_vector: FixedBitSet,
}

#[pymethods]
impl InterleavedBloomFilter {
    #[new]
    fn new(max_seq_length: usize,
           fragment_length: usize,
           overlap: usize,
           k: usize,
           num_hashes: usize,
           error_rate: f64,
           confidence: f64) -> Self {
        let num_vectors = calc_num_fragments(max_seq_length, fragment_length, overlap);
        let size_per_vector = calc_size_per_vector(num_hashes, fragment_length, k);

        let active_filter: Vec<FixedBitSet> = (0..num_vectors)
            .map(|_| {
                let mut bit_set = FixedBitSet::with_capacity(size_per_vector);
                bit_set.set_range(.., false);
                bit_set
            })
            .collect();
        let filters: HashMap<(String, String), Vec<FixedBitSet>> = HashMap::new();

        let hashes: Vec<usize> = Vec::with_capacity(num_hashes);
        let h_fn_1 = RapidInlineHasher::default();
        let h_fn_2 = RapidInlineHasher::new(0x9e3779b97f4a7c15);

        let standard_norm_dist: Normal = Normal::new(0.0, 1.0).unwrap();
        let counting_vector: Vec<usize> = vec![0; size_per_vector];
        let binning_vector = FixedBitSet::with_capacity(size_per_vector);

        Self {
            fragment_length,
            overlap,
            k,
            error_rate,
            confidence,
            num_vectors,
            size_per_vector,
            active_filter,
            filters,
            h_fn_1,
            h_fn_2,
            hashes,
            standard_norm_dist,
            counting_vector,
            binning_vector
        }
    }
    pub fn insert_sequence(&mut self, seq_id: (String, String), seq: &str) {
        if self.filters.contains_key(&seq_id) {
            return;
        }
        // NOTE: duplication
        let mut filter: Vec<FixedBitSet> = (0..self.num_vectors)
            .map(|_| {
                let mut bit_set = FixedBitSet::with_capacity(self.size_per_vector);
                bit_set.set_range(.., false);
                bit_set
            })
            .collect();

        for (i, frag_start) in (0..seq.len()).step_by(self.fragment_length - self.overlap).enumerate() {
            let fragment = &seq[frag_start..(frag_start + self.fragment_length).min(seq.len())];
            let num_k_mers = if fragment.len() >= self.k { fragment.len() - self.k } else { 0 } + 1;

            for j in 0..num_k_mers {
                hash(&mut self.h_fn_1, &mut self.h_fn_2, &mut self.hashes, fragment[j..(j + self.k).min(fragment.len())].as_bytes());
                for h in self.hashes.drain(..) {
                    FixedBitSet::set(&mut filter[h % self.num_vectors], i, true);
                }
            }
        }

        self.filters.insert(seq_id, filter);
    }

    pub fn activate_filter(&mut self, seq_id: (String, String)) {
        if ! self.filters.contains_key(&seq_id) {
            panic!("Tried to activate a filter that is not present.");
        }

        for (i, filter) in self.filters.get(&seq_id).unwrap().iter().enumerate(){
            self.active_filter[i] |= filter;
        }
    }

    pub fn is_sequence_present(&mut self, seq: &str) -> bool {
        let (_, upper) = calculate_confidence_intervals(self.standard_norm_dist, seq.len(), self.k, self.error_rate, self.confidence);
        let threshold = seq.len() - self.k + 1 - upper;
        self.is_sequence_present_aux(seq.as_bytes(), threshold) || self.is_sequence_present_aux(&reverse_complement(seq.as_bytes()), threshold)
    }
}

impl InterleavedBloomFilter {
    fn is_sequence_present_aux(&mut self, seq: &[u8], threshold: usize) -> bool {
        let mut max_count: usize = 0;
        let num_k_mers = if seq.len() >= self.k { seq.len() - self.k } else { 0 } + 1;
        self.counting_vector.fill(0);

        for i in 0..num_k_mers {
            self.binning_vector.set_range(.., true);

            hash(&mut self.h_fn_1, &mut self.h_fn_2, &mut self.hashes, &seq[i..(i + self.k).min(seq.len())]);
            for digest in self.hashes.drain(..) {
                self.binning_vector &= &self.active_filter[digest % self.num_vectors];
            }

            for j in self.binning_vector.ones() {
                self.counting_vector[j] += 1;
                max_count = max_count.max(self.counting_vector[j]);
                if max_count >= threshold {
                    return true;
                }
                if threshold - max_count > num_k_mers - i - 1 {
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
