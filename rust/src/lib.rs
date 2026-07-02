// ZOMBI2 Rust fast-path engines for the built-in UnorderedGenome + UniformRates model
// (D/T/L/O, additive uniform-recipient transfers, hard family-size cap).
//
// * simulate_profiles — genomes are per-family COPY COUNTS only (no ids); returns the
//   presence/copy-number profile matrix. Fastest; for large-scale profile studies.
// * simulate_log — tracks individual gene lineages (re-minting ids at every event, like the
//   Python engine) and emits the full event genealogy + leaf genomes, so Python can rebuild
//   the same event log, gene trees and outputs.

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use std::collections::{BTreeMap, HashMap};

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

// ===================================================================================
// Full-log engine: tracks gene lineages and emits the event genealogy (drop-in log)
// ===================================================================================

struct LGene {
    gid: u64,
    family: u64,
}

/// Order-free genome tracking individual gene lineages (plus per-family counts for O(1)
/// cap checks). Deterministic given the RNG: draws index into `genes`, never iteration order.
struct LogGenome {
    genes: Vec<LGene>,
    fam_count: HashMap<u64, u32>,
}
impl LogGenome {
    fn new() -> Self {
        LogGenome { genes: Vec::new(), fam_count: HashMap::new() }
    }
    #[inline]
    fn size(&self) -> usize {
        self.genes.len()
    }
    #[inline]
    fn copy_number(&self, family: u64) -> u32 {
        *self.fam_count.get(&family).unwrap_or(&0)
    }
    fn add(&mut self, gid: u64, family: u64) {
        self.genes.push(LGene { gid, family });
        *self.fam_count.entry(family).or_insert(0) += 1;
    }
    fn remove_at(&mut self, i: usize) -> LGene {
        let g = self.genes.swap_remove(i);
        if let Some(c) = self.fam_count.get_mut(&g.family) {
            *c -= 1;
            if *c == 0 {
                self.fam_count.remove(&g.family);
            }
        }
        g
    }
}

/// Columnar event log. Each record: event code, branch, time, donor/recipient (-1 = none),
/// family, and up to three gene ids (g1/g2 = -1 when unused). genes[0] is the incoming
/// lineage; genes[1..] the outgoing ones (the from->to genealogy edge).
struct Records {
    ev: Vec<u8>,
    br: Vec<u32>,
    tm: Vec<f64>,
    dn: Vec<i64>,
    rc: Vec<i64>,
    fm: Vec<u64>,
    g0: Vec<u64>,
    g1: Vec<i64>,
    g2: Vec<i64>,
}
impl Records {
    fn new() -> Self {
        Records { ev: vec![], br: vec![], tm: vec![], dn: vec![], rc: vec![], fm: vec![],
                  g0: vec![], g1: vec![], g2: vec![] }
    }
    #[inline]
    #[allow(clippy::too_many_arguments)]
    fn push(&mut self, ev: u8, br: u32, tm: f64, dn: i64, rc: i64, fm: u64, g0: u64, g1: i64, g2: i64) {
        self.ev.push(ev);
        self.br.push(br);
        self.tm.push(tm);
        self.dn.push(dn);
        self.rc.push(rc);
        self.fm.push(fm);
        self.g0.push(g0);
        self.g1.push(g1);
        self.g2.push(g2);
    }
}

const EV_O: u8 = 0;
const EV_D: u8 = 1;
const EV_T: u8 = 2;
const EV_L: u8 = 3;
const EV_S: u8 = 4;

struct LogEngine {
    children: Vec<Vec<usize>>,
    extant_leaf: Vec<bool>,
    d: f64,
    t: f64,
    l: f64,
    o: f64,
    cap: i64,
    rng: Rng,
    next_gid: u64,
    next_fam: u64,
    genomes: Vec<Option<LogGenome>>,
    alive_list: Vec<usize>,
    alive_pos: Vec<i64>,
    fenwick: Fenwick,
    rec: Records,
    events: u64,
    max_events: u64,
}

impl LogEngine {
    #[inline]
    fn subtotal(&self, g: &LogGenome) -> f64 {
        (g.size() as f64) * (self.d + self.t + self.l) + self.o
    }
    fn activate(&mut self, node: usize) {
        let sub = self.subtotal(self.genomes[node].as_ref().unwrap());
        self.alive_pos[node] = self.alive_list.len() as i64;
        self.alive_list.push(node);
        self.fenwick.set(node, sub);
    }
    fn refresh(&mut self, node: usize) {
        let sub = self.subtotal(self.genomes[node].as_ref().unwrap());
        self.fenwick.set(node, sub);
    }
    fn deactivate(&mut self, node: usize) -> LogGenome {
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
            self.fire(branch, t);
        }
    }

    fn fire(&mut self, branch: usize, t: f64) {
        let size = self.genomes[branch].as_ref().unwrap().size() as f64;
        let sd = size * self.d;
        let sdt = size * (self.d + self.t);
        let sdtl = size * (self.d + self.t + self.l);
        let branch_total = sdtl + self.o;
        let r = self.rng.next_f64() * branch_total;

        if r < sd {
            // duplication (skipped at the cap, after drawing the target — matches Python)
            let i = self.rng.uniform(size as u64) as usize;
            let (old, fam) = {
                let ge = &self.genomes[branch].as_ref().unwrap().genes[i];
                (ge.gid, ge.family)
            };
            if self.cap >= 0 && self.genomes[branch].as_ref().unwrap().copy_number(fam) as i64 >= self.cap {
                return;
            }
            let left = self.next_gid;
            let right = self.next_gid + 1;
            self.next_gid += 2;
            {
                let g = self.genomes[branch].as_mut().unwrap();
                g.remove_at(i);
                g.add(left, fam);
                g.add(right, fam);
            }
            self.rec.push(EV_D, branch as u32, t, -1, -1, fam, old, left as i64, right as i64);
            self.refresh(branch);
        } else if r < sdt {
            // transfer (additive; uniform recipient != donor; over-cap forces replacement losses)
            let recipient = match self.pick_recipient(branch) {
                Some(x) => x,
                None => return,
            };
            let i = self.rng.uniform(size as u64) as usize;
            let (old, fam) = {
                let ge = &self.genomes[branch].as_ref().unwrap().genes[i];
                (ge.gid, ge.family)
            };
            let cont = self.next_gid;
            let tc = self.next_gid + 1;
            self.next_gid += 2;
            {
                // donor lineage continues (re-minted in place; donor size unchanged)
                let g = self.genomes[branch].as_mut().unwrap();
                g.genes[i] = LGene { gid: cont, family: fam };
            }
            self.genomes[recipient].as_mut().unwrap().add(tc, fam);
            self.rec.push(EV_T, branch as u32, t, branch as i64, recipient as i64, fam, old, cont as i64, tc as i64);
            if self.cap >= 0 {
                let total = self.genomes[recipient].as_ref().unwrap().copy_number(fam) as i64;
                if total > self.cap {
                    for _ in 0..(total - self.cap) {
                        self.remove_one_family(recipient, fam, tc, t);
                    }
                }
            }
            self.refresh(branch);
            self.refresh(recipient);
        } else if r < sdtl {
            // loss
            let i = self.rng.uniform(size as u64) as usize;
            let gene = self.genomes[branch].as_mut().unwrap().remove_at(i);
            self.rec.push(EV_L, branch as u32, t, -1, -1, gene.family, gene.gid, -1, -1);
            self.refresh(branch);
        } else {
            // origination (brand-new family)
            let fam = self.next_fam;
            self.next_fam += 1;
            let gid = self.next_gid;
            self.next_gid += 1;
            self.genomes[branch].as_mut().unwrap().add(gid, fam);
            self.rec.push(EV_O, branch as u32, t, -1, -1, fam, gid, -1, -1);
            self.refresh(branch);
        }
    }

    /// Remove one copy of `family` from `node`, preferring a pre-existing (non-`protected`)
    /// copy; log it as a LOSS. Used for cap-forced replacement after a transfer.
    fn remove_one_family(&mut self, node: usize, family: u64, protected: u64, t: f64) {
        let pick = {
            let g = self.genomes[node].as_ref().unwrap();
            let mut pre: Vec<usize> = Vec::new();
            let mut prot: Vec<usize> = Vec::new();
            for (k, gene) in g.genes.iter().enumerate() {
                if gene.family == family {
                    if gene.gid != protected {
                        pre.push(k);
                    } else {
                        prot.push(k);
                    }
                }
            }
            let pool = if !pre.is_empty() { pre } else { prot };
            if pool.is_empty() {
                return;
            }
            pool[self.rng.uniform(pool.len() as u64) as usize]
        };
        let gene = self.genomes[node].as_mut().unwrap().remove_at(pick);
        self.rec.push(EV_L, node as u32, t, -1, -1, gene.family, gene.gid, -1, -1);
    }

    /// Speciate: re-mint each parent lineage into both children and log the bifurcation.
    fn speciate_with(&mut self, node: usize, t: f64, parent: LogGenome) {
        let c1 = self.children[node][0];
        let c2 = self.children[node][1];
        let mut g1 = LogGenome::new();
        let mut g2 = LogGenome::new();
        for gene in parent.genes.iter() {
            let n1 = self.next_gid;
            let n2 = self.next_gid + 1;
            self.next_gid += 2;
            g1.add(n1, gene.family);
            g2.add(n2, gene.family);
            self.rec.push(EV_S, node as u32, t, -1, -1, gene.family, gene.gid, n1 as i64, n2 as i64);
        }
        self.genomes[c1] = Some(g1);
        self.genomes[c2] = Some(g2);
        self.activate(c1);
        self.activate(c2);
    }
}

type LogColumns = (Vec<u8>, Vec<u32>, Vec<f64>, Vec<i64>, Vec<i64>, Vec<u64>, Vec<u64>, Vec<i64>, Vec<i64>);

#[allow(clippy::too_many_arguments)]
#[pyfunction]
fn simulate_log(
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
) -> PyResult<(LogColumns, Vec<(usize, Vec<(u64, u64)>)>)> {
    let mut children: Vec<Vec<usize>> = vec![Vec::new(); n_nodes];
    for (node, &p) in parent.iter().enumerate() {
        if p >= 0 {
            children[p as usize].push(node);
        }
    }

    let mut eng = LogEngine {
        children,
        extant_leaf,
        d: duplication,
        t: transfer,
        l: loss,
        o: origination,
        cap,
        rng: Rng::new(seed),
        next_gid: 0,
        next_fam: 0,
        genomes: (0..n_nodes).map(|_| None).collect(),
        alive_list: Vec::new(),
        alive_pos: vec![-1; n_nodes],
        fenwick: Fenwick::new(n_nodes),
        rec: Records::new(),
        events: 0,
        max_events: 2_000_000_000,
    };

    // seed the root genome (one origination per initial family), then speciate the root
    let root_t = time[root];
    let mut root_genome = LogGenome::new();
    for _ in 0..initial_size {
        let fam = eng.next_fam;
        eng.next_fam += 1;
        let gid = eng.next_gid;
        eng.next_gid += 1;
        root_genome.add(gid, fam);
        eng.rec.push(EV_O, root as u32, root_t, -1, -1, fam, gid, -1, -1);
    }
    eng.speciate_with(root, root_t, root_genome);

    let mut order: Vec<usize> = (0..n_nodes).filter(|&i| parent[i] >= 0).collect();
    order.sort_by(|&a, &b| time[a].partial_cmp(&time[b]).unwrap());

    let mut leaves: Vec<(usize, Vec<(u64, u64)>)> = Vec::new();
    let mut t = root_t;
    for node in order {
        eng.evolve_interval(t, time[node])?;
        t = time[node];
        let genome = eng.deactivate(node);
        if eng.extant_leaf[node] {
            let cols: Vec<(u64, u64)> = genome.genes.iter().map(|g| (g.gid, g.family)).collect();
            leaves.push((node, cols));
        } else if !eng.children[node].is_empty() {
            eng.speciate_with(node, t, genome);
        }
    }

    let r = eng.rec;
    Ok(((r.ev, r.br, r.tm, r.dn, r.rc, r.fm, r.g0, r.g1, r.g2), leaves))
}

#[pyfunction]
fn available() -> bool {
    true
}

#[pymodule]
fn zombi2_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(available, m)?)?;
    m.add_function(wrap_pyfunction!(simulate_profiles, m)?)?;
    m.add_function(wrap_pyfunction!(simulate_log, m)?)?;
    Ok(())
}
