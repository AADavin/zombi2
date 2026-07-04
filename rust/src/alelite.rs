//! ALElite dated-DTL reconciliation likelihood — the hot loops of `zombi2/alelite/dated.py`.
//!
//! A faithful port of the Python engine (same time-slicing, same coupled backward-Euler sweep,
//! same recursion), so the two agree to floating-point up to summation order. The Python side
//! stays the reference implementation and parity oracle; this is only the speed path.
//!
//! Everything is kept in its own module so ALElite can be lifted out of ZOMBI2 cleanly.

use pyo3::prelude::*;

const EPS: f64 = 1e-9;
const MAX_ITERS: usize = 500; // matches zombi2/alelite/undated.py
const TOL: f64 = 1e-13;

struct Slice {
    hi: f64,
    alive: Vec<usize>,
    d: usize,
    dt: f64,
    e: Vec<Vec<f64>>, // [d+1][n_alive] extinction rows, filled by compute_extinction
}

struct Engine {
    left: Vec<i64>,
    right: Vec<i64>,
    is_leaf: Vec<bool>,
    time: Vec<f64>,
    length: Vec<f64>,
    root: usize,
    d: f64,
    t: f64,
    lam: f64,
    tot: f64,
    slices: Vec<Slice>,
    e_bottom: Vec<f64>,
}

impl Engine {
    fn new(
        parent: &[i64],
        left: Vec<i64>,
        right: Vec<i64>,
        is_leaf: Vec<bool>,
        time: Vec<f64>,
        parent_time: Vec<f64>,
        root: usize,
        dup: f64,
        transfer: f64,
        loss: f64,
        n_steps: usize,
    ) -> Self {
        let n = parent.len();
        let length: Vec<f64> = (0..n).map(|i| time[i] - parent_time[i]).collect();
        let tot = dup + transfer + loss;

        // boundaries = sorted unique node times ∪ {0}
        let mut times: Vec<f64> = time.iter().map(|&x| (x * 1e12).round() / 1e12).collect();
        times.push(0.0);
        times.sort_by(|a, b| a.partial_cmp(b).unwrap());
        times.dedup();

        let mut slices = Vec::new();
        for k in 0..times.len() - 1 {
            let (lo, hi) = (times[k], times[k + 1]);
            let alive: Vec<usize> = (0..n)
                .filter(|&i| parent_time[i] <= lo + EPS && time[i] >= hi - EPS && length[i] > EPS)
                .collect();
            let span = hi - lo;
            let d = if span > 0.0 {
                n_steps.max((span * tot * 4.0).ceil() as usize)
            } else {
                1
            };
            let dt = span / (d as f64);
            slices.push(Slice { hi, alive, d, dt, e: Vec::new() });
        }

        let mut eng = Engine {
            left,
            right,
            is_leaf,
            time,
            length,
            root,
            d: dup,
            t: transfer,
            lam: loss,
            tot,
            slices,
            e_bottom: vec![0.0; n],
        };
        eng.compute_extinction();
        eng
    }

    fn compute_extinction(&mut self) {
        let (d, t, lam, tot) = (self.d, self.t, self.lam, self.tot);
        for si in (0..self.slices.len()).rev() {
            let (hi, dt, dcount, alive) = {
                let s = &self.slices[si];
                (s.hi, s.dt, s.d, s.alive.clone())
            };
            let n_al = alive.len();
            let mut e = vec![0.0f64; n_al];
            for (j, &br) in alive.iter().enumerate() {
                if (self.time[br] - hi).abs() < EPS {
                    e[j] = if self.is_leaf[br] {
                        0.0
                    } else {
                        self.e_bottom[self.left[br] as usize] * self.e_bottom[self.right[br] as usize]
                    };
                } else {
                    e[j] = self.e_bottom[br];
                }
            }
            let mut rows: Vec<Vec<f64>> = Vec::with_capacity(dcount + 1);
            rows.push(e.clone());
            for _ in 0..dcount {
                let tot_e: f64 = e.iter().sum();
                let mut new = vec![0.0; n_al];
                for j in 0..n_al {
                    let ej = e[j];
                    let ebar = if n_al > 1 { (tot_e - ej) / ((n_al - 1) as f64) } else { 0.0 };
                    let de = tot * ej - lam - d * ej * ej - t * ej * ebar;
                    new[j] = ej - dt * de;
                }
                e = new;
                rows.push(e.clone());
            }
            for (j, &br) in alive.iter().enumerate() {
                self.e_bottom[br] = e[j];
            }
            self.slices[si].e = rows;
        }
    }

    /// Log marginal likelihood of one gene tree (local post-order arrays; root = last node).
    fn gene_loglik(
        &self,
        gis_leaf: &[bool],
        gleft: &[i64],
        gright: &[i64],
        gspecies: &[i64],
        origination: u8,
    ) -> f64 {
        let (d, t, tot) = (self.d, self.t, self.tot);
        let n = self.time.len();
        let m = gis_leaf.len();
        let mut f_bottom: Vec<Vec<f64>> = (0..m).map(|_| vec![0.0; n]).collect();
        let mut f_rows: Vec<Option<Vec<Vec<Vec<f64>>>>> = (0..m).map(|_| None).collect();

        for u in 0..m {
            let internal = !gis_leaf[u];
            let (v, w) = if internal {
                (gleft[u] as usize, gright[u] as usize)
            } else {
                (usize::MAX, usize::MAX)
            };
            let mut per_slice: Vec<Vec<Vec<f64>>> = vec![Vec::new(); self.slices.len()];

            for gi in (0..self.slices.len()).rev() {
                let sl = &self.slices[gi];
                let (hi, dt, dcount) = (sl.hi, sl.dt, sl.d);
                let n_al = sl.alive.len();
                let mut f = vec![0.0f64; n_al];
                for (j, &br) in sl.alive.iter().enumerate() {
                    if (self.time[br] - hi).abs() < EPS {
                        if self.is_leaf[br] {
                            f[j] = if !internal && gspecies[u] >= 0 && gspecies[u] as usize == br {
                                1.0
                            } else {
                                0.0
                            };
                        } else {
                            let a = self.left[br] as usize;
                            let c = self.right[br] as usize;
                            let (ea, ec) = (self.e_bottom[a], self.e_bottom[c]);
                            f[j] = f_bottom[u][a] * ec + f_bottom[u][c] * ea; // SL
                            if internal {
                                f[j] += f_bottom[v][a] * f_bottom[w][c]
                                    + f_bottom[v][c] * f_bottom[w][a];
                            }
                        }
                    } else {
                        f[j] = f_bottom[u][br];
                    }
                }
                let mut rows: Vec<Vec<f64>> = Vec::with_capacity(dcount + 1);
                rows.push(f.clone());
                for s in 0..dcount {
                    let erow = &sl.e[s];
                    let tot_f: f64 = f.iter().sum();
                    let tot_e: f64 = erow.iter().sum();
                    let (fv_row, fw_row): (&[f64], &[f64]) = if internal {
                        (
                            &f_rows[v].as_ref().unwrap()[gi][s],
                            &f_rows[w].as_ref().unwrap()[gi][s],
                        )
                    } else {
                        (&[], &[])
                    };
                    let (sum_fv, sum_fw) = if internal && n_al > 1 {
                        (fv_row.iter().sum::<f64>(), fw_row.iter().sum::<f64>())
                    } else {
                        (0.0, 0.0)
                    };
                    let mut new = vec![0.0; n_al];
                    for j in 0..n_al {
                        let fj = f[j];
                        let ej = erow[j];
                        let ebar = if n_al > 1 { (tot_e - ej) / ((n_al - 1) as f64) } else { 0.0 };
                        let pbar = if n_al > 1 { (tot_f - fj) / ((n_al - 1) as f64) } else { 0.0 };
                        let homog = fj * (tot - 2.0 * d * ej - t * ebar) - t * ej * pbar;
                        let mut src = 0.0;
                        if internal {
                            let (fv, fw) = (fv_row[j], fw_row[j]);
                            src += d * fv * fw;
                            if n_al > 1 {
                                src += t / ((n_al - 1) as f64)
                                    * (fv * (sum_fw - fw) + fw * (sum_fv - fv));
                            }
                        }
                        new[j] = fj - dt * homog + dt * src;
                    }
                    f = new;
                    rows.push(f.clone());
                }
                for (j, &br) in sl.alive.iter().enumerate() {
                    f_bottom[u][br] = f[j];
                }
                per_slice[gi] = rows;
            }
            f_rows[u] = Some(per_slice);
            if internal {
                f_rows[v] = None; // children no longer needed
                f_rows[w] = None;
            }
        }

        let root_u = m - 1;
        let like = if origination == 1 {
            let reals: Vec<usize> = (0..n).filter(|&i| self.length[i] > EPS).collect();
            reals.iter().map(|&e| f_bottom[root_u][e]).sum::<f64>() / (reals.len() as f64)
        } else {
            let a = self.left[self.root] as usize;
            let c = self.right[self.root] as usize;
            let (ea, ec) = (self.e_bottom[a], self.e_bottom[c]);
            let mut like = f_bottom[root_u][a] * ec + f_bottom[root_u][c] * ea;
            if !gis_leaf[root_u] {
                let v = gleft[root_u] as usize;
                let w = gright[root_u] as usize;
                like += f_bottom[v][a] * f_bottom[w][c] + f_bottom[v][c] * f_bottom[w][a];
            }
            like
        };
        if like <= 0.0 {
            f64::NEG_INFINITY
        } else {
            like.ln()
        }
    }

    fn extinct_logprob(&self, origination: u8) -> f64 {
        let pe = if origination == 1 {
            let reals: Vec<usize> = (0..self.time.len()).filter(|&i| self.length[i] > EPS).collect();
            reals.iter().map(|&e| self.e_bottom[e]).sum::<f64>() / (reals.len() as f64)
        } else {
            self.e_bottom[self.left[self.root] as usize] * self.e_bottom[self.right[self.root] as usize]
        };
        if pe > 0.0 {
            pe.ln()
        } else {
            f64::NEG_INFINITY
        }
    }
}

// ===================================================================== undated / reldated
//
// A faithful port of zombi2/alelite/undated.py: per-branch odds (pD,pT,pL,pS), a coupled
// extinction fixed point, and a per-gene-node DP with its own SL/DL/TL fixed point. The only
// difference between undated and reldated is the transfer neighborhood `nb`: None = any branch
// (undated); an explicit per-branch list of time-overlapping branches = reldated.

struct UndatedEngine {
    n: usize,
    left: Vec<i64>,
    right: Vec<i64>,
    is_leaf: Vec<bool>,
    root: usize,
    pd: f64,
    pt: f64,
    ps: f64,
    nb: Option<Vec<Vec<usize>>>,
    e: Vec<f64>,
}

impl UndatedEngine {
    #[allow(clippy::too_many_arguments)]
    fn new(
        left: Vec<i64>,
        right: Vec<i64>,
        is_leaf: Vec<bool>,
        time: Vec<f64>,
        parent_time: Vec<f64>,
        root: usize,
        dup: f64,
        transfer: f64,
        loss: f64,
        transfers: u8,
    ) -> Self {
        let n = left.len();
        let denom = 1.0 + dup + transfer + loss;
        let (pd, pt, pl, ps) = (dup / denom, transfer / denom, loss / denom, 1.0 / denom);
        let nb = if transfers == 1 {
            let tol = 1e-9 * time.iter().cloned().fold(1.0_f64, f64::max);
            let mut nb = vec![Vec::new(); n];
            for e in 0..n {
                for f in 0..n {
                    if f == e {
                        continue;
                    }
                    let overlap = time[e].min(time[f]) - parent_time[e].max(parent_time[f]);
                    if overlap > tol {
                        nb[e].push(f);
                    }
                }
            }
            Some(nb)
        } else {
            None
        };
        let mut eng = UndatedEngine {
            n, left, right, is_leaf, root, pd, pt, ps, nb, e: vec![0.0; n],
        };
        eng.extinction(pl);
        eng
    }

    fn mean_over(&self, vec: &[f64], total: f64, e: usize) -> f64 {
        match &self.nb {
            None => {
                if self.n > 1 {
                    (total - vec[e]) / ((self.n - 1) as f64)
                } else {
                    0.0
                }
            }
            Some(nb) => {
                let lst = &nb[e];
                if lst.is_empty() {
                    0.0
                } else {
                    lst.iter().map(|&f| vec[f]).sum::<f64>() / (lst.len() as f64)
                }
            }
        }
    }

    fn extinction(&mut self, pl: f64) {
        let (pd, pt, ps, n) = (self.pd, self.pt, self.ps, self.n);
        let mut e = vec![0.0f64; n];
        for _ in 0..MAX_ITERS {
            let mut total: f64 = e.iter().sum();
            let mut delta = 0.0f64;
            for i in 0..n {
                let ebar = self.mean_over(&e, total, i);
                let child = if self.is_leaf[i] {
                    0.0
                } else {
                    e[self.left[i] as usize] * e[self.right[i] as usize]
                };
                let new = pl + pd * e[i] * e[i] + pt * e[i] * ebar + ps * child;
                delta = delta.max((new - e[i]).abs());
                total += new - e[i];
                e[i] = new;
            }
            if delta < TOL {
                break;
            }
        }
        self.e = e;
    }

    fn propagate(&self, pu: &mut [f64], a: &[f64]) {
        let (pd, pt, ps, n) = (self.pd, self.pt, self.ps, self.n);
        let total_e: f64 = self.e.iter().sum();
        for _ in 0..MAX_ITERS {
            let mut total: f64 = pu.iter().sum();
            let mut delta = 0.0f64;
            for i in 0..n {
                let ebar = self.mean_over(&self.e, total_e, i);
                let pbar = self.mean_over(pu, total, i);
                let sl = if self.is_leaf[i] {
                    0.0
                } else {
                    ps * (pu[self.left[i] as usize] * self.e[self.right[i] as usize]
                        + pu[self.right[i] as usize] * self.e[self.left[i] as usize])
                };
                let denom = 1.0 - 2.0 * pd * self.e[i] - pt * ebar;
                let new = (a[i] + sl + pt * self.e[i] * pbar) / denom;
                delta = delta.max((new - pu[i]).abs());
                total += new - pu[i];
                pu[i] = new;
            }
            if delta < TOL {
                break;
            }
        }
    }

    fn gene_loglik(&self, gis_leaf: &[bool], gleft: &[i64], gright: &[i64],
                   gspecies: &[i64], origination: u8) -> f64 {
        let (pd, pt, ps, n) = (self.pd, self.pt, self.ps, self.n);
        let m = gis_leaf.len();
        let mut p: Vec<Vec<f64>> = vec![Vec::new(); m];
        for u in 0..m {
            let mut a = vec![0.0f64; n];
            if gis_leaf[u] {
                a[gspecies[u] as usize] = ps;
            } else {
                let av = p[gleft[u] as usize].clone();
                let bv = p[gright[u] as usize].clone();
                let atot: f64 = av.iter().sum();
                let btot: f64 = bv.iter().sum();
                for e in 0..n {
                    let mut term = 2.0 * pd * av[e] * bv[e];
                    if !self.is_leaf[e] {
                        let f = self.left[e] as usize;
                        let gg = self.right[e] as usize;
                        term += ps * (av[f] * bv[gg] + av[gg] * bv[f]);
                    }
                    if pt > 0.0 {
                        let abar = self.mean_over(&av, atot, e);
                        let bbar = self.mean_over(&bv, btot, e);
                        term += pt * (av[e] * bbar + bv[e] * abar);
                    }
                    a[e] = term;
                }
            }
            let mut pu = vec![0.0f64; n];
            self.propagate(&mut pu, &a);
            p[u] = pu;
        }
        let root_u = m - 1;
        let like = if origination == 1 {
            p[root_u].iter().sum::<f64>() / (n as f64)
        } else {
            p[root_u][self.root]
        };
        if like <= 0.0 {
            f64::NEG_INFINITY
        } else {
            like.ln()
        }
    }

    fn extinct_logprob(&self, origination: u8) -> f64 {
        let pe = if origination == 1 {
            self.e.iter().sum::<f64>() / (self.n as f64)
        } else {
            self.e[self.root]
        };
        if pe > 0.0 {
            pe.ln()
        } else {
            f64::NEG_INFINITY
        }
    }
}

/// Joint undated/reldated log-likelihood of a batch of gene trees. `transfers` is `0`
/// (undated — any recipient) or `1` (reldated — only time-overlapping recipients).
#[pyfunction]
#[allow(clippy::too_many_arguments)]
pub fn undated_joint_loglik(
    sp_left: Vec<i64>,
    sp_right: Vec<i64>,
    sp_is_leaf: Vec<bool>,
    sp_time: Vec<f64>,
    sp_parent_time: Vec<f64>,
    sp_root: usize,
    gt_offsets: Vec<usize>,
    gt_is_leaf: Vec<bool>,
    gt_left: Vec<i64>,
    gt_right: Vec<i64>,
    gt_species: Vec<i64>,
    dup: f64,
    transfer: f64,
    loss: f64,
    transfers: u8,
    origination: u8,
    n_extinct: usize,
) -> f64 {
    let eng = UndatedEngine::new(
        sp_left, sp_right, sp_is_leaf, sp_time, sp_parent_time, sp_root,
        dup, transfer, loss, transfers,
    );
    let mut ll = 0.0f64;
    if n_extinct > 0 {
        ll += (n_extinct as f64) * eng.extinct_logprob(origination);
    }
    for i in 0..gt_offsets.len().saturating_sub(1) {
        let (a, b) = (gt_offsets[i], gt_offsets[i + 1]);
        ll += eng.gene_loglik(
            &gt_is_leaf[a..b], &gt_left[a..b], &gt_right[a..b], &gt_species[a..b], origination,
        );
    }
    ll
}

/// Joint dated log-likelihood of a batch of gene trees sharing one species tree and rates.
///
/// The species tree is post-order flat arrays; gene trees are concatenated (CSR `gt_offsets`),
/// each with local post-order `is_leaf`/`left`/`right` (child = local index, `-1` for a leaf)
/// and `species` (the species-tree leaf-branch index for a tip, `-1` internal). `origination`
/// is `0` (root) or `1` (uniform); `n_extinct` adds `k·log P(no survivor)`.
#[pyfunction]
#[allow(clippy::too_many_arguments)]
pub fn dated_joint_loglik(
    sp_parent: Vec<i64>,
    sp_left: Vec<i64>,
    sp_right: Vec<i64>,
    sp_is_leaf: Vec<bool>,
    sp_time: Vec<f64>,
    sp_parent_time: Vec<f64>,
    sp_root: usize,
    gt_offsets: Vec<usize>,
    gt_is_leaf: Vec<bool>,
    gt_left: Vec<i64>,
    gt_right: Vec<i64>,
    gt_species: Vec<i64>,
    dup: f64,
    transfer: f64,
    loss: f64,
    n_steps: usize,
    origination: u8,
    n_extinct: usize,
) -> f64 {
    let eng = Engine::new(
        &sp_parent, sp_left, sp_right, sp_is_leaf, sp_time, sp_parent_time, sp_root,
        dup, transfer, loss, n_steps,
    );
    let mut ll = 0.0f64;
    if n_extinct > 0 {
        ll += (n_extinct as f64) * eng.extinct_logprob(origination);
    }
    for i in 0..gt_offsets.len().saturating_sub(1) {
        let (a, b) = (gt_offsets[i], gt_offsets[i + 1]);
        ll += eng.gene_loglik(
            &gt_is_leaf[a..b],
            &gt_left[a..b],
            &gt_right[a..b],
            &gt_species[a..b],
            origination,
        );
    }
    ll
}

/// Per-family dated log-likelihoods: extinction is built once and every gene tree in the batch
/// is scored against it, returning one log-lik per tree (in input order) — for a per-family
/// score table rather than the joint sum.
#[pyfunction]
#[allow(clippy::too_many_arguments)]
pub fn dated_family_loglik(
    sp_parent: Vec<i64>,
    sp_left: Vec<i64>,
    sp_right: Vec<i64>,
    sp_is_leaf: Vec<bool>,
    sp_time: Vec<f64>,
    sp_parent_time: Vec<f64>,
    sp_root: usize,
    gt_offsets: Vec<usize>,
    gt_is_leaf: Vec<bool>,
    gt_left: Vec<i64>,
    gt_right: Vec<i64>,
    gt_species: Vec<i64>,
    dup: f64,
    transfer: f64,
    loss: f64,
    n_steps: usize,
    origination: u8,
) -> Vec<f64> {
    let eng = Engine::new(
        &sp_parent, sp_left, sp_right, sp_is_leaf, sp_time, sp_parent_time, sp_root,
        dup, transfer, loss, n_steps,
    );
    (0..gt_offsets.len().saturating_sub(1))
        .map(|i| {
            let (a, b) = (gt_offsets[i], gt_offsets[i + 1]);
            eng.gene_loglik(&gt_is_leaf[a..b], &gt_left[a..b], &gt_right[a..b], &gt_species[a..b], origination)
        })
        .collect()
}

/// Per-family undated/reldated log-likelihoods (extinction built once). `transfers` = 0
/// (undated) or 1 (reldated).
#[pyfunction]
#[allow(clippy::too_many_arguments)]
pub fn undated_family_loglik(
    sp_left: Vec<i64>,
    sp_right: Vec<i64>,
    sp_is_leaf: Vec<bool>,
    sp_time: Vec<f64>,
    sp_parent_time: Vec<f64>,
    sp_root: usize,
    gt_offsets: Vec<usize>,
    gt_is_leaf: Vec<bool>,
    gt_left: Vec<i64>,
    gt_right: Vec<i64>,
    gt_species: Vec<i64>,
    dup: f64,
    transfer: f64,
    loss: f64,
    transfers: u8,
    origination: u8,
) -> Vec<f64> {
    let eng = UndatedEngine::new(
        sp_left, sp_right, sp_is_leaf, sp_time, sp_parent_time, sp_root, dup, transfer, loss, transfers,
    );
    (0..gt_offsets.len().saturating_sub(1))
        .map(|i| {
            let (a, b) = (gt_offsets[i], gt_offsets[i + 1]);
            eng.gene_loglik(&gt_is_leaf[a..b], &gt_left[a..b], &gt_right[a..b], &gt_species[a..b], origination)
        })
        .collect()
}
