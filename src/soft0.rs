// soft scattering: 
//   rustc soft0.rs
// examples:
//   ./soft0
//   ./soft0 sphere-0.01905-66-128.msh
// MM 14.4.2026

use std::f64::consts::PI;

// ---------- Complex arithmetic ----------

type C = (f64, f64); // (re, im)

fn cadd(x: C, y: C) -> C { (x.0 + y.0, x.1 + y.1) }
fn csub(x: C, y: C) -> C { (x.0 - y.0, x.1 - y.1) }
fn cmul(x: C, y: C) -> C { (x.0 * y.0 - x.1 * y.1, x.0 * y.1 + x.1 * y.0,)}
fn cdiv(x: C, y: C) -> C { let m = y.0 * y.0 + y.1 * y.1;
    ( (x.0 * y.0 + x.1 * y.1) / m, (x.1 * y.0 - x.0 * y.1) / m,)
}
fn cabs(z: C) -> f64 {(z.0 * z.0 + z.1 * z.1).sqrt()}
fn cexpj(x: f64) -> C {(x.cos(), x.sin())}

// ---------- Gauss elimination (complex) ----------

fn gauss(a: &[C], b: &[C]) -> Vec<C> {
    let n = b.len();
    let mut mat = a.to_vec();
    let mut vec = b.to_vec();   
    for i in 0..n {                   // 1. Forward elimination
        let pivot = mat[i * n + i];
        for j in i..n {               // normalize row i
            mat[i * n + j] = cdiv(mat[i * n + j], pivot);
        }
        vec[i] = cdiv(vec[i], pivot);
        for k in (i + 1)..n {         // clear below pivot
            let factor = mat[k * n + i];
            for j in i..n {
                let subtrahend = cmul(factor, mat[i * n + j]);
                mat[k * n + j] = csub(mat[k * n + j], subtrahend);
            }
            let v_subtrahend = cmul(factor, vec[i]);
            vec[k] = csub(vec[k], v_subtrahend);
        }
    }    
    let mut x = vec![ (0.0, 0.0); n ];// 2. Back substitution
    for i in (0..n).rev() {
        let mut sum = (0.0, 0.0);
        for j in (i + 1)..n {
            sum = cadd(sum, cmul(mat[i * n + j], x[j]));
        }
        x[i] = csub(vec[i], sum);
    }
    x
}

// ---------- Geometry / mesh ----------

#[derive(Clone, Debug)]
struct MeshData {
    x: Vec<f64>, y: Vec<f64>, z: Vec<f64>,
    facets_0: Vec<usize>, facets_1: Vec<usize>, facets_2: Vec<usize>,
}

fn centroids(d: &MeshData) -> Vec<(f64, f64, f64)> {
    let n = d.facets_0.len();
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let i1 = d.facets_0[i];
        let i2 = d.facets_1[i];
        let i3 = d.facets_2[i];
        let xc = d.x[i1] + d.x[i2] + d.x[i3];
        let yc = d.y[i1] + d.y[i2] + d.y[i3];
        let zc = d.z[i1] + d.z[i2] + d.z[i3];
        out.push((-xc / 3.0, -zc / 3.0, -yc / 3.0));  // rotated version!!!
    }
    out
}

fn dist(x: (f64, f64, f64), y: (f64, f64, f64)) -> f64 {
    let dx = x.0 - y.0; let dy = x.1 - y.1; let dz = x.2 - y.2;
    (dx * dx + dy * dy + dz * dz).sqrt()
}
fn dot(x: (f64, f64, f64), y: (f64, f64, f64)) -> f64 { x.0 * y.0 + x.1 * y.1 + x.2 * y.2}

// ---------- Green’s functions ----------

fn green(x: &[(f64, f64, f64)], y: &[(f64, f64, f64)], k: f64) -> Vec<C> {
    let n = x.len();
    let mut g = Vec::with_capacity(n * n);
    for i in 0..n {
        for j in 0..n {
            let r = dist(x[i], y[j]);
            if r == 0.0 {
                g.push((0.0, k / (4.0 * PI)));
            } else {
                g.push(cdiv(cexpj(k * r), (4.0 * PI * r, 0.0)));
            }
        }
    }
    g
}

fn green2(x: &[(f64, f64, f64)], y: &[(f64, f64, f64)], k: f64) -> Vec<C> {
    let mut out = Vec::with_capacity(y.len() * x.len());
    for v in y {
        for xj in x {
            let phase = -k * dot(*xj, *v);
            out.push((phase.cos() / (4.0 * PI), phase.sin() / (4.0 * PI)));
        }
    }
    out
}

// ---------- Matrix–vector product (complex) ----------

fn matvec(a: &[C], x: &[C]) -> Vec<C> {
    let n = x.len();
    let m = a.len() / n;
    let mut y = Vec::with_capacity(m);
    for i in 0..m {
        let mut sum = (0.0, 0.0);
        for j in 0..n {
            sum = cadd(sum, cmul(a[i * n + j], x[j]));
        }
        y.push(sum);
    }
    y
}

// ---------- Directions, incident field ----------

fn dir(th: f64, ph: f64) -> (f64, f64, f64) {
    let th0 = th.to_radians() + PI;
    let ph0 = ph.to_radians() + PI;
    let st = th0.sin();
    let ct = th0.cos();
    let sp = ph0.sin();
    let cp = ph0.cos();
    (ct * cp, st, -ct * sp)
}

fn pinc(x: &[(f64, f64, f64)], k: f64, d: (f64, f64, f64)) -> Vec<C> {
    let mut b = Vec::with_capacity(x.len());
    for xi in x {
        let phase = k * dot(*xi, d);
        b.push(cmul(cexpj(phase), (-1.0, 0.0)));
    }
    b
}

// ---------- Fixed formatting (like Typst fixed) ----------

/*
fn fixed(num: f64, digits: usize) -> String {
    let s = format!("{:.1$}", num, digits);
    s
}
*/

// --- sphere ----

use std::collections::HashMap;

fn midpoint(i: usize, j: usize, r: f64, vx: &mut Vec<[f64; 3]>, cache: &mut HashMap<(usize, usize), usize>, ) -> usize {
    let key = if i < j { (i, j) } else { (j, i) };
    if let Some(&idx) = cache.get(&key) {
        return idx;
    }
    let p = [ (vx[i][0]+vx[j][0])*0.5, (vx[i][1]+vx[j][1])*0.5, (vx[i][2]+vx[j][2])*0.5,];
    let l = (p[0]*p[0] + p[1]*p[1] + p[2]*p[2]).sqrt();
    let p = [p[0]/l * r, p[1]/l * r, p[2]/l * r];
    let idx = vx.len();
    vx.push(p);
    cache.insert(key, idx);
    idx
}

pub fn sphere(r: f64, n: usize) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let mut vx: Vec<[f64; 3]> = vec![         // initial octahedron
        [0.0, 0.0,  r], [0.0, 0.0, -r], [-r, 0.0, 0.0],
        [0.0, -r, 0.0], [ r, 0.0, 0.0], [0.0,  r, 0.0],
    ];
    let mut fx: Vec<[usize; 3]> = vec![
        [0, 3, 4], [0, 4, 5], [0, 5, 2], [0, 2, 3],
        [1, 4, 3], [1, 5, 4], [1, 2, 5], [1, 3, 2],
    ];
    let mut midpoint_cache: HashMap<(usize, usize), usize> = HashMap::new();
    for _ in 0..n {     // subdivision loop
        let mut new_fx = Vec::new();
        midpoint_cache.clear();
        for &[i1, i2, i3] in &fx {
            let a = midpoint(i1, i2, r, &mut vx, &mut midpoint_cache);
            let b = midpoint(i2, i3, r, &mut vx, &mut midpoint_cache);
            let c = midpoint(i3, i1, r, &mut vx, &mut midpoint_cache);
            new_fx.push([i1, a, c]);
            new_fx.push([a, i2, b]);
            new_fx.push([c, b, i3]);
            new_fx.push([a, b, c]);
        }
        fx = new_fx;
    }
    (vx, fx)
}
/*
//---
fn norm(v: [f64; 3]) -> [f64; 3] {
    let l = (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]).sqrt();
    [v[0]/l, v[1]/l, v[2]/l]
}

fn sphere(r: f64, n: usize) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let mut vx = vec![
        [0.0, 0.0,  r], [0.0, 0.0, -r], [-r, 0.0, 0.0],
        [0.0, -r, 0.0], [ r, 0.0, 0.0], [0.0,  r, 0.0],
    ];
    let mut fx = vec![
        [0, 3, 4],[0, 4, 5],[0, 5, 2],[0, 2, 3],[1, 4, 3],[1, 5, 4],[1, 2, 5],[1, 3, 2],
    ];
    let mut nv = vx.len();
    let mut nt = fx.len();
    for _ in 0..n {
        let ntold = nt;
        for i in 0..ntold {
            let [i1, i2, i3] = fx[i];
            let a0 = [
                (vx[i1][0] + vx[i2][0]) * 0.5,
                (vx[i1][1] + vx[i2][1]) * 0.5,
                (vx[i1][2] + vx[i2][2]) * 0.5,
            ];
            let b0 = [
                (vx[i2][0] + vx[i3][0]) * 0.5,
                (vx[i2][1] + vx[i3][1]) * 0.5,
                (vx[i2][2] + vx[i3][2]) * 0.5,
            ];
            let c0 = [
                (vx[i3][0] + vx[i1][0]) * 0.5,
                (vx[i3][1] + vx[i1][1]) * 0.5,
                (vx[i3][2] + vx[i1][2]) * 0.5,
            ];
            let a = norm(a0);
            let b = norm(b0);
            let c = norm(c0);
            vx.push([r * a[0], r * a[1], r * a[2]]);
            vx.push([r * b[0], r * b[1], r * b[2]]);
            vx.push([r * c[0], r * c[1], r * c[2]]);
            let a_i = nv;
            let b_i = nv + 1;
            let c_i = nv + 2;
            fx.push([i1, a_i, c_i]);
            fx.push([a_i, i2, b_i]);
            fx.push([b_i, i3, c_i]);
            fx[i] = [a_i, b_i, c_i];
            nv += 3;
            nt += 3;
        }
    }
    (vx, fx)
}
*/

//---

fn mesh_from_vx_fx(vx: &Vec<[f64; 3]>, fx: &Vec<[usize; 3]>) -> MeshData {
    let mut x = Vec::with_capacity(vx.len());
    let mut y = Vec::with_capacity(vx.len());
    let mut z = Vec::with_capacity(vx.len());
    for v in vx {
        x.push(v[0]); y.push(v[1]); z.push(v[2]);
    }
    let mut facets_0 = Vec::with_capacity(fx.len());
    let mut facets_1 = Vec::with_capacity(fx.len());
    let mut facets_2 = Vec::with_capacity(fx.len());
    for f in fx {
        facets_0.push(f[0]); facets_1.push(f[1]); facets_2.push(f[2]);
    }
    MeshData {x, y, z, facets_0, facets_1, facets_2,}
}

// ---- file reader

use std::fs;

pub fn msh_from_file(path: &str) -> std::io::Result<(Vec<[f64; 3]>, Vec<[usize; 3]>)> {
    let text = fs::read_to_string(path)?;
    Ok(msh(&text))
}

pub fn msh(s: &str) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let lines: Vec<&str> = s.lines().collect();
    let n: usize = lines[4].trim().parse().unwrap();     // number of vertices
    let mut v = Vec::with_capacity(n);                   // parse vertices
    for line in &lines[5 .. 5 + n] {
        let parts: Vec<&str> = line.trim().split_whitespace().collect();
        let x: f64 = parts[1].parse().unwrap();
        let y: f64 = parts[2].parse().unwrap();
        let z: f64 = parts[3].parse().unwrap();
        v.push([x, y, z]);
    }
    let m: usize = lines[n + 7].trim().parse().unwrap(); // number of faces
    let mut e = Vec::with_capacity(m);   // parse faces (convert to 0-based indexing)
    for line in &lines[n + 8 .. n + 8 + m] {
        let parts: Vec<&str> = line.trim().split_whitespace().collect();
        let a: usize = parts[5].parse::<usize>().unwrap() - 1;
        let b: usize = parts[6].parse::<usize>().unwrap() - 1;
        let c: usize = parts[7].parse::<usize>().unwrap() - 1;
        e.push([a, b, c]);
    }
    (v, e)
}


// ---------- main ----------
use std::env;
use std::time::Instant;

fn main() {
    let args: Vec<String> = env::args().collect();
    let (vx, fx);
    let (mut f0, mut c0, mut th, mut ph) = (38e3_f64, 1480.0_f64, 270.0_f64, 0.0_f64);
    if args.len() > 1 {
	eprintln!("{}:",args[1]);
	(vx, fx) = msh_from_file(&args[1]).unwrap(); // "sphere-0.01905-66-128.msh"
    } else {
	let r = 0.01905;
        eprintln!("r={}:",r);
        (vx, fx) = sphere(r, 2);
    }
    if args.len() > 2 { f0 = args[2].parse().unwrap();}
    if args.len() > 3 { c0 = args[3].parse().unwrap();}
    if args.len() > 4 { th = args[4].parse().unwrap();}
    if args.len() > 5 { ph = args[5].parse().unwrap();}

/*
    for (i, v) in vx.iter().enumerate() {
      println!("{:4}: [{:.6}, {:.6}, {:.6}]", i, v[0], v[1], v[2]);
    }
    for (i, f) in fx.iter().enumerate() {
      println!("{:4}: [{}, {}, {}]", i, f[0], f[1], f[2]);
    }
*/
    let mesh = mesh_from_vx_fx(&vx, &fx);

    let x = centroids(&mesh);
    let k = 2.0 * PI * f0 / c0;

    let a = green(&x, &x, k);
    let d = dir(th, ph);
    let b = pinc(&x, k, d);
    let start = Instant::now();
    let q = gauss(&a, &b);
    let elapsed = start.elapsed(); 

    let mut y_dirs = Vec::with_capacity(360);
    for t in 0..360 {
        y_dirs.push(dir(-(t as f64 + 180.0), ph));
    }

    let s = green2(&x, &y_dirs, k);
    let psc = matvec(&s, &q);

    let idx = ((th + 180.0) % 360.0) as usize;
    let lbs = cabs(psc[idx]);

    for i in 0..psc.len() {
      println!("{} {:.6}", i, cabs(psc[i]));
    }

    eprintln!("lbs = {}", lbs);
    eprintln!("Execution time: {:.3}", elapsed.as_secs_f64());
}
