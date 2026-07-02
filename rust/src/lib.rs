// ZOMBI2 Rust fast-path engine (increment 1): profiles-only, independent gene families.
//
// A genome is just per-family COPY COUNTS (no gene ids / trees / event log yet), so this
// mirrors the Python UnorderedGenome + UniformRates model (D/T/L/O, additive uniform
// transfers, hard family-size cap) and returns the presence/copy-number profile matrix.
// The full event log / gene trees remain the Python engine.

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use std::collections::BTreeMap;

/// xoshiro256** — small, fast, deterministic PRNG (seeded via splitmix64).
struct Rng {
    s: [u64; 4],
}
impl Rng {
    fn new(seed: u64) -> Self {
        let mut z = seed;
        let mut sm = || {
            z = z.wrapping_add(0x9E3779B97F4A7C15);
            let mut x = z;
            x = (x ^ (x >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94D049BB133111EB);
            x ^ (x >> 31)
        };
        Rng { s: [sm(), sm(), sm(), sm()] }
    }
    #[inline]
    fn next_u64(&mut self) -> u64 {
        let result = self.s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        result
    }
    #[inline]
    fn next_f64(&mut self) -> f64 {
        // 53-bit mantissa in [0, 1)
        (self.next_u64() >> 11) as f64 * (1.0 / 9007199254740992.0)
    }
    #[inline]
    fn uniform(&mut self, n: u64) -> u64 {
        self.next_u64() % n
    }
    #[inline]
    fn exponential(&mut self, rate: f64) -> f64 {
        // 1 - u lies in (0, 1], so ln is safe
        -(1.0 - self.next_f64()).ln() / rate
    }
}

/// Fenwick tree over branch subtotals (0-indexed slots; value 0 means "not a live branch").
struct Fenwick {
    n: usize,
    tree: Vec<f64>,
    vals: Vec<f64>,
    total: f64,
}
impl Fenwick {
    fn new(n: usize) -> Self {
        Fenwick { n, tree: vec![0.0; n + 1], vals: vec![0.0; n], total: 0.0 }
    }
    fn set(&mut self, i: usize, v: f64) {
        let d = v - self.vals[i];
        if d == 0.0 {
            return;
        }
        self.vals[i] = v;
        self.total += d;
        let mut j = i + 1;
        while j <= self.n {
            self.tree[j] += d;
            j += j & j.wrapping_neg();
        }
    }
    fn find(&self, mut value: f64) -> usize {
        let mut pos = 0usize;
        let mut k = 1usize;
        while (k << 1) <= self.n {
            k <<= 1;
        }
        while k > 0 {
            let nxt = pos + k;
            if nxt <= self.n && self.tree[nxt] < value {
                pos = nxt;
                value -= self.tree[nxt];
            }
            k >>= 1;
        }
        pos
    }
}

/// Order-free genome: family id -> copy count (BTreeMap for deterministic iteration).
#[derive(Clone)]
struct Genome {
    counts: BTreeMap<u64, u32>,
    size: u64,
}
impl Genome {
    fn new() -> Self {
        Genome { counts: BTreeMap::new(), size: 0 }
    }
    fn add(&mut self, fam: u64) {
        *self.counts.entry(fam).or_insert(0) += 1;
        self.size += 1;
    }
    fn remove_one(&mut self, fam: u64) {
        if let Some(c) = self.counts.get_mut(&fam) {
            *c -= 1;
            self.size -= 1;
            if *c == 0 {
                self.counts.remove(&fam);
            }
        }
    }
    fn copy_number(&self, fam: u64) -> u32 {
        *self.counts.get(&fam).unwrap_or(&0)
    }
    /// The family holding the r-th gene copy (0 <= r < size), in key order.
    fn family_at(&self, mut r: u64) -> u64 {
        for (&fam, &c) in self.counts.iter() {
            let c = c as u64;
            if r < c {
                return fam;
            }
            r -= c;
        }
        unreachable!("family_at index out of range")
    }
}

struct Engine {
    children: Vec<Vec<usize>>,
    extant_leaf: Vec<bool>,
    d: f64,
    t: f64,
    l: f64,
    o: f64,
    cap: i64, // -1 = no cap
    rng: Rng,
    next_family: u64,
    genomes: Vec<Option<Genome>>,
    alive_list: Vec<usize>,
    alive_pos: Vec<i64>, // node -> index in alive_list, -1 if not alive
    fenwick: Fenwick,
    events: u64,
    max_events: u64,
}

impl Engine {
    fn subtotal(&self, g: &Genome) -> f64 {
        (g.size as f64) * (self.d + self.t + self.l) + self.o
    }

    fn activate(&mut self, node: usize, genome: Genome) {
        let sub = self.subtotal(&genome);
        self.genomes[node] = Some(genome);
        self.alive_pos[node] = self.alive_list.len() as i64;
        self.alive_list.push(node);
        self.fenwick.set(node, sub);
    }

    fn deactivate(&mut self, node: usize) -> Genome {
        self.fenwick.set(node, 0.0);
        let pos = self.alive_pos[node] as usize;
        let last = *self.alive_list.last().unwrap();
        self.alive_list.swap_remove(pos);
        if last != node {
            self.alive_pos[last] = pos as i64;
        }
        self.alive_pos[node] = -1;
        self.genomes[node].take().unwrap()
    }

    fn refresh(&mut self, node: usize) {
        let sub = self.subtotal(self.genomes[node].as_ref().unwrap());
        self.fenwick.set(node, sub);
    }

    /// Pick a uniform recipient among alive branches other than `donor` (None if <2 alive).
    fn pick_recipient(&mut self, donor: usize) -> Option<usize> {
        let k = self.alive_list.len();
        if k < 2 {
            return None;
        }
        let mut j = self.rng.uniform((k - 1) as u64) as usize;
        let donor_pos = self.alive_pos[donor] as usize;
        if j >= donor_pos {
            j += 1;
        }
        Some(self.alive_list[j])
    }

    fn evolve_interval(&mut self, mut t: f64, t1: f64) -> PyResult<()> {
        loop {
            let total = self.fenwick.total;
            if total <= 0.0 {
                return Ok(());
            }
            let dt = self.rng.exponential(total);
            if !dt.is_finite() || t + dt >= t1 {
                return Ok(());
            }
            t += dt;
            self.events += 1;
            if self.events > self.max_events {
                return Err(PyRuntimeError::new_err(
                    "exceeded max_events; families likely growing without bound — set max_family_size",
                ));
            }
            let branch = self.fenwick.find((1.0 - self.rng.next_f64()) * total);
            self.fire(branch);
        }
    }

    fn fire(&mut self, branch: usize) {
        let (size, sd, sdt, sdtl) = {
            let g = self.genomes[branch].as_ref().unwrap();
            let size = g.size as f64;
            (size, size * self.d, size * (self.d + self.t), size * (self.d + self.t + self.l))
        };
        let branch_total = sdtl + self.o;
        let r = self.rng.next_f64() * branch_total;

        if r < sd {
            // duplication
            let g = self.genomes[branch].as_mut().unwrap();
            let idx = self.rng.uniform(size as u64);
            let fam = g.family_at(idx);
            if self.cap < 0 || (g.copy_number(fam) as i64) < self.cap {
                g.add(fam);
            }
            self.refresh(branch);
        } else if r < sdt {
            // transfer (additive; over-cap recipient is a net-zero no-op)
            let fam = {
                let g = self.genomes[branch].as_ref().unwrap();
                let idx = self.rng.uniform(size as u64);
                g.family_at(idx)
            };
            if let Some(recipient) = self.pick_recipient(branch) {
                let rg = self.genomes[recipient].as_mut().unwrap();
                if self.cap < 0 || (rg.copy_number(fam) as i64) < self.cap {
                    rg.add(fam);
                    self.refresh(recipient);
                }
            }
        } else if r < sdtl {
            // loss
            let g = self.genomes[branch].as_mut().unwrap();
            let idx = self.rng.uniform(size as u64);
            let fam = g.family_at(idx);
            g.remove_one(fam);
            self.refresh(branch);
        } else {
            // origination
            let fam = self.next_family;
            self.next_family += 1;
            self.genomes[branch].as_mut().unwrap().add(fam);
            self.refresh(branch);
        }
    }
}

#[allow(clippy::too_many_arguments)]
#[pyfunction]
fn simulate_profiles(
    n_nodes: usize,
    parent: Vec<i64>,
    time: Vec<f64>,
    extant_leaf: Vec<bool>,
    root: usize,
    duplication: f64,
    transfer: f64,
    loss: f64,
    origination: f64,
    initial_size: usize,
    cap: i64,
    seed: u64,
) -> PyResult<Vec<(usize, Vec<(u64, u32)>)>> {
    // children adjacency from parent pointers
    let mut children: Vec<Vec<usize>> = vec![Vec::new(); n_nodes];
    for (node, &p) in parent.iter().enumerate() {
        if p >= 0 {
            children[p as usize].push(node);
        }
    }

    let mut engine = Engine {
        children,
        extant_leaf,
        d: duplication,
        t: transfer,
        l: loss,
        o: origination,
        cap,
        rng: Rng::new(seed),
        next_family: 0,
        genomes: (0..n_nodes).map(|_| None).collect(),
        alive_list: Vec::new(),
        alive_pos: vec![-1; n_nodes],
        fenwick: Fenwick::new(n_nodes),
        events: 0,
        max_events: 2_000_000_000,
    };

    // seed the root genome, then speciate the root into its children
    let mut root_genome = Genome::new();
    for _ in 0..initial_size {
        let fam = engine.next_family;
        engine.next_family += 1;
        root_genome.add(fam);
    }
    for &child in engine.children[root].clone().iter() {
        engine.activate(child, root_genome.clone());
    }

    // node events (all non-root nodes) in time order
    let mut order: Vec<usize> = (0..n_nodes).filter(|&i| parent[i] >= 0).collect();
    order.sort_by(|&a, &b| time[a].partial_cmp(&time[b]).unwrap());

    let mut profiles: Vec<(usize, Vec<(u64, u32)>)> = Vec::new();
    let mut t = time[root];
    for node in order {
        engine.evolve_interval(t, time[node])?;
        t = time[node];
        let genome = engine.deactivate(node);
        if engine.extant_leaf[node] {
            let cols: Vec<(u64, u32)> = genome.counts.iter().map(|(&f, &c)| (f, c)).collect();
            profiles.push((node, cols));
        } else {
            for &child in engine.children[node].clone().iter() {
                engine.activate(child, genome.clone());
            }
        }
    }
    Ok(profiles)
}

#[pyfunction]
fn available() -> bool {
    true
}

#[pymodule]
fn zombi2_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(available, m)?)?;
    m.add_function(wrap_pyfunction!(simulate_profiles, m)?)?;
    Ok(())
}
