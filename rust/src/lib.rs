// ZOMBI2 Rust fast-path engines for the built-in UnorderedGenome + UniformRates model
// (D/T/L/O, additive uniform-recipient transfers, hard family-size cap).
//
// * simulate_profiles — genomes are per-family COPY COUNTS only (no ids); returns the
//   presence/copy-number profile matrix. Fastest; for large-scale profile studies.
// * simulate_log — tracks individual gene lineages (re-minting ids at every event, like the
//   Python engine) and emits the full event genealogy + leaf genomes, so Python can rebuild
//   the same event log, gene trees and outputs.

use pyo3::exceptions::{PyIOError, PyRuntimeError};
use pyo3::prelude::*;
use rayon::prelude::*;
use std::collections::{BTreeMap, HashMap};
use std::fmt::Write as _;
use std::fs::{self, File};
use std::io::{BufWriter, Write as _};
use std::path::Path;

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

/// Transfer mechanics (mirrors Python's TransferModel). `decay < 0` means uniform recipient.
#[derive(Clone, Copy)]
struct Transfers {
    replacement: f64,
    decay: f64,
    allow_self: bool,
}

/// Choose a transfer recipient among the alive branches: uniform, or weighted by
/// exp(-decay * (d - dmin)) with patristic distance d = 2*(t - t_MRCA). Optionally allows
/// the donor itself (a self-transfer = duplication). Returns None if no candidate exists.
#[allow(clippy::too_many_arguments)]
fn choose_recipient(
    rng: &mut Rng,
    alive_list: &[usize],
    alive_pos: &[i64],
    donor: usize,
    parent: &[i64],
    time: &[f64],
    t: f64,
    tr: &Transfers,
) -> Option<usize> {
    let k = alive_list.len();
    if k == 0 || (!tr.allow_self && k < 2) {
        return None;
    }
    if tr.decay < 0.0 {
        // uniform recipient
        if tr.allow_self {
            return Some(alive_list[rng.uniform(k as u64) as usize]);
        }
        let donor_pos = alive_pos[donor] as usize;
        let mut j = rng.uniform((k - 1) as u64) as usize;
        if j >= donor_pos {
            j += 1;
        }
        return Some(alive_list[j]);
    }
    // distance-weighted: mark the donor's ancestor chain, then walk each candidate up to it
    let mut anc: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let mut node = donor as i64;
    while node >= 0 {
        anc.insert(node as usize);
        node = parent[node as usize];
    }
    let mut dists: Vec<f64> = Vec::with_capacity(k);
    for &r in alive_list {
        if r == donor {
            dists.push(if tr.allow_self { 0.0 } else { f64::INFINITY });
            continue;
        }
        let mut n = r;
        while !anc.contains(&n) {
            n = parent[n] as usize;
        }
        dists.push(2.0 * (t - time[n]));
    }
    let dmin = dists.iter().cloned().fold(f64::INFINITY, f64::min);
    if !dmin.is_finite() {
        return None;
    }
    let mut total = 0.0;
    let weights: Vec<f64> = dists
        .iter()
        .map(|&d| {
            let w = if d.is_finite() { (-tr.decay * (d - dmin)).exp() } else { 0.0 };
            total += w;
            w
        })
        .collect();
    if total <= 0.0 {
        return None;
    }
    let mut rr = rng.next_f64() * total;
    for (idx, &w) in weights.iter().enumerate() {
        rr -= w;
        if rr <= 0.0 {
            return Some(alive_list[idx]);
        }
    }
    Some(*alive_list.last().unwrap())
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
    parent: Vec<i64>,
    time: Vec<f64>,
    tr: Transfers,
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
            // transfer: additive; recipient by TransferModel; then cap-forced or
            // probabilistic replacement (drop a pre-existing copy). Self-transfer = duplication.
            let fam = {
                let g = self.genomes[branch].as_ref().unwrap();
                let idx = self.rng.uniform(size as u64);
                g.family_at(idx)
            };
            let recipient = match choose_recipient(&mut self.rng, &self.alive_list,
                &self.alive_pos, branch, &self.parent, &self.time, t, &self.tr) {
                Some(x) => x,
                None => return,
            };
            self.genomes[recipient].as_mut().unwrap().add(fam);
            let total = self.genomes[recipient].as_ref().unwrap().copy_number(fam) as i64;
            let mut removals = 0i64;
            if self.cap >= 0 && total > self.cap {
                removals = total - self.cap;
            } else if total - 1 >= 1 && self.tr.replacement > 0.0 && self.rng.next_f64() < self.tr.replacement {
                removals = 1;
            }
            for _ in 0..removals {
                self.genomes[recipient].as_mut().unwrap().remove_one(fam);
            }
            self.refresh(recipient);
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
    replacement: f64,
    distance_decay: f64,
    allow_self: bool,
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
        parent: parent.clone(),
        time: time.clone(),
        tr: Transfers { replacement, decay: distance_decay, allow_self },
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
    parent: Vec<i64>,
    time: Vec<f64>,
    tr: Transfers,
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
            // transfer: recipient by TransferModel (uniform / distance-weighted / self);
            // donor lineage bifurcates (old -> cont in donor, tc into recipient); then
            // cap-forced or probabilistic replacement drops a pre-existing recipient copy.
            let recipient = match choose_recipient(&mut self.rng, &self.alive_list,
                &self.alive_pos, branch, &self.parent, &self.time, t, &self.tr) {
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
            let total = self.genomes[recipient].as_ref().unwrap().copy_number(fam) as i64;
            let mut removals = 0i64;
            if self.cap >= 0 && total > self.cap {
                removals = total - self.cap;
            } else if total - 1 >= 1 && self.tr.replacement > 0.0 && self.rng.next_f64() < self.tr.replacement {
                removals = 1;
            }
            for _ in 0..removals {
                self.remove_one_family(recipient, fam, tc, t);
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

/// Run the full-log engine; return its records and extant leaf genomes.
#[allow(clippy::too_many_arguments)]
fn run_log(
    n_nodes: usize,
    parent: &[i64],
    time: &[f64],
    extant_leaf: &[bool],
    root: usize,
    duplication: f64,
    transfer: f64,
    loss: f64,
    origination: f64,
    initial_size: usize,
    cap: i64,
    seed: u64,
    tr: Transfers,
) -> PyResult<(Records, Vec<(usize, Vec<(u64, u64)>)>)> {
    let mut children: Vec<Vec<usize>> = vec![Vec::new(); n_nodes];
    for (node, &p) in parent.iter().enumerate() {
        if p >= 0 {
            children[p as usize].push(node);
        }
    }

    let mut eng = LogEngine {
        children,
        extant_leaf: extant_leaf.to_vec(),
        parent: parent.to_vec(),
        time: time.to_vec(),
        tr,
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

    Ok((eng.rec, leaves))
}

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
    replacement: f64,
    distance_decay: f64,
    allow_self: bool,
) -> PyResult<(LogColumns, Vec<(usize, Vec<(u64, u64)>)>)> {
    let tr = Transfers { replacement, decay: distance_decay, allow_self };
    let (r, leaves) = run_log(n_nodes, &parent, &time, &extant_leaf, root,
        duplication, transfer, loss, origination, initial_size, cap, seed, tr)?;
    Ok(((r.ev, r.br, r.tm, r.dn, r.rc, r.fm, r.g0, r.g1, r.g2), leaves))
}

// ============ output writing in Rust (gene trees + tables) — no big return to Python ============

/// %g-like formatting: `sig` significant figures with trailing zeros trimmed.
fn fmt_g(x: f64, sig: i32) -> String {
    if x == 0.0 || !x.is_finite() {
        return "0".to_string();
    }
    let exp = x.abs().log10().floor() as i32;
    let decimals = (sig - 1 - exp).clamp(0, 17) as usize;
    let s = format!("{:.*}", decimals, x);
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        s
    }
}

/// Natural-ish sort key matching Python's _natkey: (numeric run in the name, name).
fn natkey(name: &str) -> (i128, &str) {
    let digits: String = name.chars().filter(|c| c.is_ascii_digit()).collect();
    (digits.parse::<i128>().unwrap_or(0), name)
}

const EV_CHAR: [&str; 5] = ["O", "D", "T", "L", "S"];
const ROLES: [&[&str]; 5] = [
    &["origin"],
    &["parent", "left", "right"],
    &["parent", "donor_copy", "transfer_copy"],
    &["lost"],
    &["parent", "child", "child"],
];

/// A reconstructed gene-tree node (species = leaf node index for extant tips).
struct GNode {
    gid: u64,
    birth: f64,
    end: f64,
    children: Vec<GNode>,
    is_loss: bool,
    species: Option<usize>,
}

fn bl(n: &GNode) -> f64 {
    (n.end - n.birth).max(0.0)
}

fn to_newick(n: &GNode, names: &[String], out: &mut String) {
    if n.children.is_empty() {
        if n.is_loss {
            let _ = write!(out, "LOSS_{}:{}", n.gid, fmt_g(bl(n), 6));
        } else if let Some(li) = n.species {
            let _ = write!(out, "{}_{}:{}", names[li], n.gid, fmt_g(bl(n), 6));
        } else {
            let _ = write!(out, "{}:{}", n.gid, fmt_g(bl(n), 6));
        }
    } else {
        out.push('(');
        for (i, c) in n.children.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            to_newick(c, names, out);
        }
        let _ = write!(out, "){}:{}", n.gid, fmt_g(bl(n), 6));
    }
}

/// Prune to lineages leading to an extant leaf; suppress degree-two nodes.
fn prune(mut n: GNode) -> Option<GNode> {
    if n.children.is_empty() {
        return if n.species.is_some() { Some(n) } else { None };
    }
    let mut kept: Vec<GNode> = n.children.into_iter().filter_map(prune).collect();
    match kept.len() {
        0 => None,
        1 => {
            let mut survivor = kept.pop().unwrap();
            survivor.birth = n.birth;
            Some(survivor)
        }
        _ => {
            n.children = kept;
            Some(n)
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn build_gnode(
    gid: u64,
    children: &HashMap<u64, (u64, u64)>,
    losses: &std::collections::HashSet<u64>,
    end_time: &HashMap<u64, f64>,
    birth: &HashMap<u64, f64>,
    gid2leaf: &HashMap<u64, usize>,
    total_age: f64,
) -> GNode {
    let b = *birth.get(&gid).unwrap_or(&0.0);
    if let Some(&(c1, c2)) = children.get(&gid) {
        GNode {
            gid, birth: b, end: end_time[&gid], is_loss: false, species: None,
            children: vec![
                build_gnode(c1, children, losses, end_time, birth, gid2leaf, total_age),
                build_gnode(c2, children, losses, end_time, birth, gid2leaf, total_age),
            ],
        }
    } else if losses.contains(&gid) {
        GNode { gid, birth: b, end: end_time[&gid], children: vec![], is_loss: true, species: None }
    } else {
        GNode { gid, birth: b, end: total_age, children: vec![], is_loss: false, species: gid2leaf.get(&gid).copied() }
    }
}

/// Reconstruct (complete, extant) Newick for one family from its time-ordered records.
fn reconstruct(
    idxs: &[usize],
    rec: &Records,
    gid2leaf: &HashMap<u64, usize>,
    names: &[String],
    total_age: f64,
) -> (Option<String>, Option<String>) {
    let mut children: HashMap<u64, (u64, u64)> = HashMap::new();
    let mut losses: std::collections::HashSet<u64> = std::collections::HashSet::new();
    let mut end_time: HashMap<u64, f64> = HashMap::new();
    let mut birth: HashMap<u64, f64> = HashMap::new();
    let mut root: Option<u64> = None;

    for &ri in idxs {
        let ev = rec.ev[ri];
        let tm = rec.tm[ri];
        match ev {
            EV_O => {
                root = Some(rec.g0[ri]);
                birth.entry(rec.g0[ri]).or_insert(tm);
            }
            EV_D | EV_T | EV_S => {
                let frm = rec.g0[ri];
                let c1 = rec.g1[ri] as u64;
                let c2 = rec.g2[ri] as u64;
                children.insert(frm, (c1, c2));
                end_time.insert(frm, tm);
                birth.insert(c1, tm);
                birth.insert(c2, tm);
            }
            EV_L => {
                losses.insert(rec.g0[ri]);
                end_time.insert(rec.g0[ri], tm);
            }
            _ => {}
        }
    }

    let root = match root {
        Some(r) => r,
        None => return (None, None),
    };
    let node = build_gnode(root, &children, &losses, &end_time, &birth, gid2leaf, total_age);
    let mut complete = String::new();
    to_newick(&node, names, &mut complete);
    complete.push(';');
    let extant = prune(node).map(|n| {
        let mut s = String::new();
        to_newick(&n, names, &mut s);
        s.push(';');
        s
    });
    (Some(complete), extant)
}

/// Write all gene-family outputs; return (n_families, n_species).
fn write_outputs(
    rec: &Records,
    leaves: &[(usize, Vec<(u64, u64)>)],
    names: &[String],
    total_age: f64,
    outdir: &str,
) -> std::io::Result<(usize, usize)> {
    let base = Path::new(outdir);
    fs::create_dir_all(base.join("gene_family_events"))?;
    fs::create_dir_all(base.join("gene_trees"))?;
    let n_events = rec.ev.len();

    // gid -> leaf index, and per-leaf family counts
    let mut gid2leaf: HashMap<u64, usize> = HashMap::new();
    let mut leaf_counts: HashMap<usize, HashMap<u64, u32>> = HashMap::new();
    for (li, pairs) in leaves {
        let cmap = leaf_counts.entry(*li).or_default();
        for &(gid, fam) in pairs {
            gid2leaf.insert(gid, *li);
            *cmap.entry(fam).or_insert(0) += 1;
        }
    }

    // species columns in natural-name order
    let mut species_cols: Vec<usize> = leaves.iter().map(|(li, _)| *li).collect();
    species_cols.sort_by(|&a, &b| natkey(&names[a]).cmp(&natkey(&names[b])));
    let n_species = species_cols.len();

    // group record indices by family (records already time-ordered)
    let mut fam_records: HashMap<u64, Vec<usize>> = HashMap::new();
    for i in 0..n_events {
        fam_records.entry(rec.fm[i]).or_default().push(i);
    }
    let mut families: Vec<u64> = fam_records.keys().copied().collect();
    families.sort_unstable();
    let n_families = families.len();

    // per-family extant copy totals (family -> (total copies, #species))
    let mut fam_copies: HashMap<u64, (u64, u64)> = HashMap::new();
    for m in leaf_counts.values() {
        for (&fam, &c) in m {
            let e = fam_copies.entry(fam).or_insert((0, 0));
            e.0 += c as u64;
            e.1 += 1;
        }
    }

    // per-family: event table + gene trees. Families are independent, so reconstruct and
    // write them in parallel (each thread writes its own files, reads shared data).
    families.par_iter().try_for_each(|&fam| -> std::io::Result<()> {
        let idxs = &fam_records[&fam];
        let label = fam + 1;

        let mut ev_s = String::from("time\tevent\tbranch\tdonor\trecipient\tnodes\n");
        for &ri in idxs {
            let code = rec.ev[ri] as usize;
            let roles = ROLES[code];
            let mut nodes = format!("{}={}", roles[0], rec.g0[ri]);
            if roles.len() == 3 {
                let _ = write!(nodes, ";{}={};{}={}", roles[1], rec.g1[ri], roles[2], rec.g2[ri]);
            }
            let donor = if rec.dn[ri] >= 0 { names[rec.dn[ri] as usize].as_str() } else { "" };
            let recip = if rec.rc[ri] >= 0 { names[rec.rc[ri] as usize].as_str() } else { "" };
            let _ = writeln!(ev_s, "{}\t{}\t{}\t{}\t{}\t{}",
                fmt_g(rec.tm[ri], 10), EV_CHAR[code], names[rec.br[ri] as usize], donor, recip, nodes);
        }
        fs::write(base.join("gene_family_events").join(format!("{}_events.tsv", label)), ev_s)?;

        let (complete, extant) = reconstruct(idxs, rec, &gid2leaf, names, total_age);
        if let Some(c) = complete {
            fs::write(base.join("gene_trees").join(format!("{}_complete.nwk", label)), c + "\n")?;
        }
        if let Some(e) = extant {
            fs::write(base.join("gene_trees").join(format!("{}_extant.nwk", label)), e + "\n")?;
        }
        Ok(())
    })?;

    // Transfers.tsv
    let mut tr = String::from("time\tfamily\tdonor_branch\trecipient_branch\tparent_id\tdonor_copy_id\ttransfer_id\n");
    for i in 0..n_events {
        if rec.ev[i] == EV_T {
            let _ = writeln!(tr, "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                fmt_g(rec.tm[i], 10), rec.fm[i] + 1,
                names[rec.dn[i] as usize], names[rec.rc[i] as usize],
                rec.g0[i], rec.g1[i], rec.g2[i]);
        }
    }
    fs::write(base.join("Transfers.tsv"), tr)?;

    // Gene_family_summary.tsv
    let mut sm = String::from("family\torigin_time\torigin_branch\tn_dup\tn_transfer\tn_loss\tn_speciation\textant_copies\tspecies_present\n");
    for &fam in &families {
        let (mut nd, mut nt, mut nl, mut ns) = (0u64, 0u64, 0u64, 0u64);
        let mut ot = String::new();
        let mut ob = String::new();
        for &ri in &fam_records[&fam] {
            match rec.ev[ri] {
                EV_D => nd += 1,
                EV_T => nt += 1,
                EV_L => nl += 1,
                EV_S => ns += 1,
                EV_O => {
                    ot = fmt_g(rec.tm[ri], 10);
                    ob = names[rec.br[ri] as usize].clone();
                }
                _ => {}
            }
        }
        let (copies, sp) = fam_copies.get(&fam).copied().unwrap_or((0, 0));
        let _ = writeln!(sm, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            fam + 1, ot, ob, nd, nt, nl, ns, copies, sp);
    }
    fs::write(base.join("Gene_family_summary.tsv"), sm)?;

    // Profiles.tsv + Presence.tsv (families x species), streamed (large at scale)
    let mut ext_fams: Vec<u64> = fam_copies.keys().copied().collect();
    ext_fams.sort_unstable();
    let col_maps: Vec<&HashMap<u64, u32>> = species_cols.iter().map(|li| &leaf_counts[li]).collect();

    let mut header = String::from("family");
    for &li in &species_cols {
        header.push('\t');
        header.push_str(&names[li]);
    }
    header.push('\n');

    let mut pf = BufWriter::new(File::create(base.join("Profiles.tsv"))?);
    let mut pr = BufWriter::new(File::create(base.join("Presence.tsv"))?);
    pf.write_all(header.as_bytes())?;
    pr.write_all(header.as_bytes())?;
    // Format rows in parallel (each family is independent), in bounded batches to cap memory,
    // then write each batch in order.
    for batch in ext_fams.chunks(512) {
        let rows: Vec<(String, String)> = batch
            .par_iter()
            .map(|&fam| {
                let mut prow = String::new();
                let mut rrow = String::new();
                let _ = write!(prow, "{}", fam + 1);
                let _ = write!(rrow, "{}", fam + 1);
                for m in &col_maps {
                    let c = *m.get(&fam).unwrap_or(&0);
                    let _ = write!(prow, "\t{}", c);
                    rrow.push('\t');
                    rrow.push(if c > 0 { '1' } else { '0' });
                }
                prow.push('\n');
                rrow.push('\n');
                (prow, rrow)
            })
            .collect();
        for (prow, rrow) in &rows {
            pf.write_all(prow.as_bytes())?;
            pr.write_all(rrow.as_bytes())?;
        }
    }
    pf.flush()?;
    pr.flush()?;

    Ok((n_families, n_species))
}

#[allow(clippy::too_many_arguments)]
#[pyfunction]
fn simulate_and_write(
    n_nodes: usize,
    parent: Vec<i64>,
    time: Vec<f64>,
    extant_leaf: Vec<bool>,
    root: usize,
    names: Vec<String>,
    duplication: f64,
    transfer: f64,
    loss: f64,
    origination: f64,
    initial_size: usize,
    cap: i64,
    seed: u64,
    replacement: f64,
    distance_decay: f64,
    allow_self: bool,
    total_age: f64,
    outdir: String,
) -> PyResult<(usize, usize, usize)> {
    let tr = Transfers { replacement, decay: distance_decay, allow_self };
    let (rec, leaves) = run_log(n_nodes, &parent, &time, &extant_leaf, root,
        duplication, transfer, loss, origination, initial_size, cap, seed, tr)?;
    let (n_families, n_species) = write_outputs(&rec, &leaves, &names, total_age, &outdir)
        .map_err(|e| PyIOError::new_err(e.to_string()))?;
    Ok((n_families, rec.ev.len(), n_species))
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
    m.add_function(wrap_pyfunction!(simulate_and_write, m)?)?;
    Ok(())
}
