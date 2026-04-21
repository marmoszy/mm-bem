// mm-bem functions for typst
// MM 15.4.2026

#import calc: max, min, round, rem, sqrt, pow, sin, cos, pi

// ---- variables
#let W = 10cm
#let zm = 32
#let (a, b, c) =  ( 0.4, -0.2, 0 )
#let a11 = -cos(c)*sin(a) - sin(c)*cos(a)*sin(b)
#let a12 =  cos(c)*cos(a) - sin(c)*sin(a)*sin(b)
#let a13 = -sin(c)*cos(b)
#let a21 =  cos(c)*cos(a)*sin(b) - sin(c)*sin(a)
#let a22 =  cos(c)*sin(a)*sin(b) + sin(c)*cos(a)
#let a23 =  cos(c)*cos(b)
#let P = (x, y, z) => (a11*z+a12*x+a13*y, a21*z+a22*x+a23*y) 
#let V = (x, y) => ((1+zm*x)*W/2, (1-zm*y)*W/2)

// ---- code ---------------

#let cadd(x,y) = (x.at(0)+y.at(0), x.at(1)+y.at(1))
#let csub(x,y) = (x.at(0)-y.at(0), x.at(1)-y.at(1))
#let cmul(x, y) = (
  x.at(0) * y.at(0) - x.at(1) * y.at(1),
  x.at(0) * y.at(1) + x.at(1) * y.at(0)
)
#let cdiv(x, y) = {
  let m = y.at(0) * y.at(0) + y.at(1) * y.at(1)
  (
    (x.at(0) * y.at(0) + x.at(1) * y.at(1)) / m,
    (x.at(1) * y.at(0) - x.at(0) * y.at(1)) / m
  )
}
#let cabs(z) = {
  sqrt(pow(z.at(0), 2) + pow(z.at(1), 2))
}
#let cexpj(x) = (cos(x),sin(x))
#let gauss(A, b) = {
  let n = b.len()
  let mat = A   // copy 
  let vec = b   // copy 
  for i in range(n) { // 1. Forward subst - triangulation
    let pivot = mat.at(i * n + i) // diagonal
    for j in range(i, n) { // Normalization of i row
      mat.at(i * n + j) = cdiv(mat.at(i * n + j), pivot)
    }
    vec.at(i) = cdiv(vec.at(i), pivot)
    for k in range(i + 1, n) { // Clear col under pivot
      let factor = mat.at(k * n + i)
      for j in range(i, n) {
        let subtrahend = cmul(factor, mat.at(i * n + j))
        mat.at(k * n + j) = csub(mat.at(k * n + j), subtrahend)
      }
      let v_subtrahend = cmul(factor, vec.at(i))
      vec.at(k) = csub(vec.at(k), v_subtrahend)
    }
  }  
  let x = range(n).map(_ => (0.0, 0.0)) 
  for i in range(n - 1, -1, step: -1) { // 2. Backsubst
    let sum = (0.0, 0.0)
    for j in range(i + 1, n) {
      sum = cadd(sum, cmul(mat.at(i * n + j), x.at(j)))
    }
    x.at(i) = csub(vec.at(i), sum)
  }
  return x
}
#let centroids(d) = {
  let m = d.x.len()
  let n = d.facets_0.len()
  let x = ()
  for i in range(n) {
    let i1 = d.facets_0.at(i)
    let i2 = d.facets_1.at(i)
    let i3 = d.facets_2.at(i)
    let xc = d.x.at(i1) + d.x.at(i2) + d.x.at(i3)
    let yc = d.y.at(i1) + d.y.at(i2) + d.y.at(i3)
    let zc = d.z.at(i1) + d.z.at(i2) + d.z.at(i3)
    //x.push((xc/3, yc/3, zc/3))
    x.push((-xc/3, -zc/3, -yc/3)) // rotate
  }
  return x
}
#let dist(x,y) = {
  sqrt(pow(x.at(0)-y.at(0),2)+pow(x.at(1)-y.at(1),2)+pow(x.at(2)-y.at(2),2))
}
#let dot(x,y) = {
  	x.at(0)*y.at(0) + x.at(1)*y.at(1) + x.at(2)*y.at(2)
}
#let green(x, y, k) = {
  let G = ()
  let n = x.len()
  for i in range(n) {
    for j in range(n) {
      let r = dist(x.at(i), y.at(j))
      if r==0 {
        G.push((0, k/(4*pi))) 
      } else {
        G.push(cdiv(cexpj(k*r), (4*pi*r, 0)))
      }
    }
  }
  return G
}
#let green2(x, y, k) = {
  let pi = calc.pi
  y.map(v=>x.map(xj=>
    (calc.cos(-k*dot(xj,v))/(4*pi), calc.sin(-k*dot(xj,v))/(4*pi)))
  ).reduce((a, b) => a + b)
}
#let matvec(A, x) = {
  let m = int(A.len() / x.len())
  let n = x.len()
  let y = ()
  for i in range(m) {
    let sum = (0.0, 0.0)
    for j in range(n) {
      sum = cadd(sum, cmul(A.at(i * n + j), x.at(j))) // sum += a_ij * x_j
    }
    y.push(sum)
  } 
  return y
}
#let dir(th, ph) = {
  let th0 = th*pi/180+pi
  let ph0 = ph*pi/180+pi
  let st = sin(th0)
  let ct = cos(th0)
  let sp = sin(ph0)
  let cp = cos(ph0)
  (ct*cp, st, -ct*sp)
}

#let pinc(x, k, d) = {
  let b = ()
  let n = x.len()
  for i in range(n) {
    b.push(cmul(cexpj(k*dot(x.at(i), d)), (-1,0)))
  }
  return b
}
#let fixed(num, digits) = {
  let s = str(round(num, digits: digits))
  let parts = s.split(".")
  if parts.len() == 1 {
    s + "." + "0" * digits
  } else {
    s + "0" * (digits - parts.at(1).len())
  }
}

// ---- sphere -----

#let midpoint(p1, p2, r) = {
  let m = (
    (p1.at(0) + p2.at(0)) / 2,
    (p1.at(1) + p2.at(1)) / 2,
    (p1.at(2) + p2.at(2)) / 2,
  )
  let len = calc.sqrt(
    m.at(0)*m.at(0) + m.at(1)*m.at(1) + m.at(2)*m.at(2)
  )
  let s = r / len
  (m.at(0)*s, m.at(1)*s, m.at(2)*s)
}
#let norm(p) = {
  let l = sqrt(p.at(0) * p.at(0) + p.at(1) * p.at(1) + p.at(2) * p.at(2))
  (p.at(0) / l, p.at(1) / l, p.at(2) / l)
}

#let sphere(r, n) = {
  let vx = (
    (0, 0,  r),(0, 0, -r),(-r, 0, 0),(0, -r, 0),( r, 0, 0),(0,  r, 0),
  ) 
  let fx = (
    (0, 3, 4),(0, 4, 5),(0, 5, 2),(0, 2, 3),
    (1, 4, 3),(1, 5, 4),(1, 2, 5),(1, 3, 2),
  )
/* icosphere
  let phi = (1 + sqrt(5)) / 2
  let vx = ( 
    (-1,  phi,  0),  ( 1,  phi,  0),
    (-1, -phi,  0),  ( 1, -phi,  0),
    ( 0, -1,  phi),  ( 0,  1,  phi), 
    ( 0, -1, -phi),  ( 0,  1, -phi),
    ( phi, 0, -1),   ( phi, 0, 1),  
    (-phi, 0, -1),   (-phi, 0, 1) 
  )
  let fx = ( 
    (0, 11, 5),   (0, 5, 1),  
    (0, 1, 7),    (0, 7, 10), 
    (0, 10, 11),  (1, 5, 9), 
    (5, 11, 4),   (11, 10, 2),
    (10, 7, 6),   (7, 1, 8),
    (3, 9, 4),    (3, 4, 2),  
    (3, 2, 6),    (3, 6, 8), 
    (3, 8, 9),    (4, 9, 5), 
    (2, 4, 11),   (6, 2, 10), 
    (8, 6, 7),    (9, 8, 1),
  )
*/
  for _ in range(n) {
    let new_fx = ()
    let cache = (:)
    for f in fx {
      let i1 = f.at(0)
      let i2 = f.at(1)
      let i3 = f.at(2)
      // edge (i1, i2)
      let key12 = str(min(i1, i2)) + "-" + str(max(i1, i2))
      let a = if cache.keys().contains(key12) {
        cache.at(key12)
      } else {
        vx.push(midpoint(vx.at(i1), vx.at(i2), r) ) 
        let idx = vx.len() - 1
        cache.insert(key12, idx)
        idx
      }
      // edge (i2, i3)
      let key23 = str(calc.min(i2, i3)) + "-" + str(calc.max(i2, i3))
      let b = if cache.keys().contains(key23) {
        cache.at(key23)
      } else {
        vx.push(midpoint(vx.at(i2), vx.at(i3), r))
        let idx = vx.len() - 1
        cache.insert(key23, idx)
        idx
      }
      // edge (i3, i1)
      let key31 = str(calc.min(i3, i1)) + "-" + str(calc.max(i3, i1))
      let c = if cache.keys().contains(key31) {
        cache.at(key31)
      } else {
        vx.push(midpoint(vx.at(i3), vx.at(i1), r))
        let idx = vx.len() - 1
        cache.insert(key31, idx)
        idx
      }
      new_fx.push((i1, a, c))
      new_fx.push((a, i2, b))
      new_fx.push((c, b, i3))
      new_fx.push((a, b, c))
    }
    fx = new_fx
  }
  (vx, fx)
}

// ----- sphere0

#let norm(v) = {
  let l = sqrt(v.at(0)*v.at(0) + v.at(1)*v.at(1) + v.at(2)*v.at(2))
  (v.at(0)/l, v.at(1)/l, v.at(2)/l)
}

#let sphere0(r, n) = {
  let vx = (
    (0, 0,  r),(0, 0, -r),(-r, 0, 0),(0, -r, 0),( r, 0, 0),(0,  r, 0),
  ) 
  let fx = (
    (0, 3, 4),(0, 4, 5),(0, 5, 2),(0, 2, 3),
    (1, 4, 3),(1, 5, 4),(1, 2, 5),(1, 3, 2),
  )
  let nv = vx.len()
  let nt = fx.len()
  for _ in range(n) {
    let ntold = nt
    for i in range(ntold) {
      let f = fx.at(i)
      let i1 = f.at(0)
      let i2 = f.at(1)
      let i3 = f.at(2)
      let a0 = (
        (vx.at(i1).at(0) + vx.at(i2).at(0)) / 2,
        (vx.at(i1).at(1) + vx.at(i2).at(1)) / 2,
        (vx.at(i1).at(2) + vx.at(i2).at(2)) / 2,
      )
      let b0 = (
        (vx.at(i2).at(0) + vx.at(i3).at(0)) / 2,
        (vx.at(i2).at(1) + vx.at(i3).at(1)) / 2,
        (vx.at(i2).at(2) + vx.at(i3).at(2)) / 2,
      )
      let c0 = (
        (vx.at(i3).at(0) + vx.at(i1).at(0)) / 2,
        (vx.at(i3).at(1) + vx.at(i1).at(1)) / 2,
        (vx.at(i3).at(2) + vx.at(i1).at(2)) / 2,
      )
      let a = norm(a0)
      let b = norm(b0)
      let c = norm(c0)
      vx.push((r*a.at(0), r*a.at(1), r*a.at(2)))
      vx.push((r*b.at(0), r*b.at(1), r*b.at(2)))
      vx.push((r*c.at(0), r*c.at(1), r*c.at(2)))
      fx.push((i1,  nv, nv+2))
      fx.push((nv,  i2, nv+1))
      fx.push((nv+1,i3, nv+2))
      fx.at(i) = (nv, nv+1, nv+2)
      nv = nv + 3
      nt = nt + 3
    }
  }
  (vx, fx)
}

// ---- vx,fx from file

#let vx-fx-from-json(fname,i,j) = {
  let j = json(fname).specimens.at(i).shapes.at(j)
  (j.x.zip(j.y,j.z), j.facets_0.zip(j.facets_1,j.facets_2))
}

#let vx-fx-from-msh(fname) = {
  let msh = csv(fname)
  let (vx, fx)  = ( (), () ) 
  let n = int(msh.at(4).at(0))
  for i in range(n) {
    let arr = msh.at(5+i).at(0).split().map(p=>float(p))
    vx.push(arr.slice(1,4))
  }
  let m = int(msh.at(4+n+3).at(0))
  for i in range(m) {
    let arr = msh.at(5+n+3+i).at(0).split().map(p=>int(p)-1)
    fx.push(arr.slice(5,8))
  }
  (vx, fx)
}

#let mesh-from-vx-fx(vx, fx) = {
  // allocate coordinate arrays
  let x = ()
  let y = ()
  let z = ()
  for v in vx {
    x.push(v.at(0))
    y.push(v.at(1))
    z.push(v.at(2))
  }
  let facets_0 = ()
  let facets_1 = ()
  let facets_2 = ()
  for f in fx {
    facets_0.push(f.at(0))
    facets_1.push(f.at(1))
    facets_2.push(f.at(2))
  }
  ( x: x, y: y, z: z, facets_0: facets_0, facets_1: facets_1, facets_2: facets_2,)
}

#let bbox(vx) = {
  let (minx,maxx) = (vx.at(0).at(0),vx.at(0).at(0))
  let (miny,maxy) = (vx.at(0).at(1),vx.at(0).at(1))
  let (minz,maxz) = (vx.at(0).at(2),vx.at(0).at(2))
  for i in range(vx.len()) {
    if(minx > vx.at(i).at(0)) { minx = vx.at(i).at(0)}
    if(maxx < vx.at(i).at(0)) { maxx = vx.at(i).at(0)}
    if(miny > vx.at(i).at(1)) { miny = vx.at(i).at(1)}
    if(maxy < vx.at(i).at(1)) { maxy = vx.at(i).at(1)}
    if(minz > vx.at(i).at(2)) { minz = vx.at(i).at(2)}
    if(maxz < vx.at(i).at(2)) { maxz = vx.at(i).at(2)}
  }
  (minx,maxx,miny,maxy,minz,maxz)
}

// ---- 3D conversion

#let mesh-points(vx,fx) = {
  let pts = ()
  for i in range(fx.len()) {
    let v1 = V(..P(..vx.at(fx.at(i).at(0))))
    let v2 = V(..P(..vx.at(fx.at(i).at(1))))
    let v3 = V(..P(..vx.at(fx.at(i).at(2))))
    pts.push( (v1,v2,v3) )
  }
  pts
}

#let line-points(vx) = {
  let pts = ()
  for i in range(vx.len()) {
    let v = V(..P(..vx.at(i)))
    pts.push(v)
  }
  pts
}

#let box-points(minx,maxx,miny,maxy,minz,maxz) = {
  let l1 = ()
  for m in (minx,2*minx/3+maxx/3,minx/3+2*maxx/3,maxx) {
    l1.push(line-points( 
      ( (m,miny,minz),(m,miny,maxz),
        (m,maxy,maxz),(m,maxy,minz) ),
    ))
  }
  for m in (miny,2*miny/3+maxy/3,miny/3+2*maxy/3,maxy) {
    l1.push(line-points( 
      ( (minx,m,minz),(minx,m,maxz),
        (maxx,m,maxz),(maxx,m,minz) ),
    ))
  }
  for m in (minz,2*minz/3+maxz/3,minz/3+2*maxz/3,maxz) {
    l1.push(line-points( 
      ( (minx,miny,m),(minx,maxy,m),
        (maxx,maxy,m),(maxx,miny,m) ),
    ))
  }
  l1
}

// ---- drawing 

#let place-points(pts,stroke) = {
  for i in range(pts.len()) {
    place(left + top, 
      curve(
        stroke: stroke,
        curve.move(pts.at(i).at(0)),
        ..pts.at(i).slice(1).map(p => curve.line(p)),
        curve.close()
      )
    )   
  }
}
