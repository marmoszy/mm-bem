#set heading(numbering: "1.1.")
#set math.equation(numbering: "(1)")
#set par(justify: true)
#set text(size: 12pt)
#show figure.where(kind: table): set figure.caption(position: top)

#title("mm-bem2")

This is typst example of using bem code inside typesetting document. Here below you will find code for visualizing generated mesh treated later as input data for BEM algorithm and output polar plot and textual representation of scattering length function.

#set text(size: 12pt)


#import "mm-code.typ": *
#let r = 0.01905
#let res = 1
#let (vx, fx) = sphere(r, res)

= Input: sphere(#r,#res)

// ----- processing 

#let x = centroids(mesh-from-vx-fx(vx, fx))
#let (f0,c0,th,ph) = (38e3,1480,270,0)
#let k = 2*calc.pi*f0/c0
#let A = green(x, x, k)
#let d = dir(th, ph)
#let b = pinc(x, k, d)
#let Q = gauss(A, b)
#let y = range(360).map(t=>dir(-(t+180), ph))
#let S = green2(x, y, k)
#let psc = matvec(S, Q)
#let lbs = cabs(psc.at(calc.rem(th+180,360)))

// ---- drawing mesh

\
\
#let ax = line-points( ((0,0,0),(2*r,0,0)) )
#let ay = line-points( ((0,0,0),(0,2*r,0)) )
#let az = line-points( ((0,0,0),(0,0,4*r)) )
#let mp = mesh-points(vx, fx)
#let bp = box-points(..bbox(vx))

#align(center)[
  #box(width: W, height: W, {
    place-points((ax,),color.red)
    place-points((ay,),color.green)
    place-points((az,),color.blue)  
    place-points(mp,color.blue)
    place-points(bp,color.black+0.2pt)
  })
]



= Output: polar plot

\
#set text(size: 12pt)

#let Pts = ()
#let W = 15cm
#let off = 0.5cm
#let pts = ()
#let sc = W/2/50
#for i in range(psc.len()) {
  let y = cabs(psc.at(i))
  let y2 = 20 * calc.log(y) + 60
  if (y2 < 0) { y2 = 0 } 
  let px = (sc * y2 + off) * calc.cos(-i*1deg)
  let py = (sc * y2 + off) * calc.sin(-i*1deg)
   pts.push((px + W/2, py + W/2))
}
#Pts.push(pts)


// ---- polar plot
#align(center)[
 #box(width: W, height: W, {
  for i in range(0,5) {
    place(center + horizon, circle(radius: off+i*W/10, stroke: gray))    
  }
  for angle in range(0, 180, step: 30) {
    place(center + horizon, 
      rotate(angle * 1deg, 
        line(
          start: (0cm, 0cm), 
          end: (W/5*(5-1) + 2*off, 0cm), 
          stroke: 0.5pt + color.gray
        )
      ) 
    )
  } 
  for angle in range(0, 360, step: 30) {
    let st = angle * 1deg
    let r2 = W/5/2*(5-1)+off+0.4cm
    let x = calc.cos(-st) * r2 
    let y = calc.sin(-st) * r2    
    place(center + horizon, dx: x, dy: y, [
      #align(center + horizon)[#angle°]
    ])
  }
  for i in range(5) {
    place(center+horizon,
      dx:off+W/10*i, dy:0.3cm, [#(-60+10*i)])
    place(center+horizon,
      dx:-(off+W/10*i), dy:0.3cm, [#(-60+10*i)])
  }
  let colors = ( color.red, color.green, color.blue)
  for (i, pts) in Pts.enumerate() {
    let color = colors.at(rem(i,colors.len()))
    place(left + top, 
    //path(..psc, stroke: red + 1.5pt)
      curve(
        stroke: color + 1.5pt,
        fill: color.transparentize(95%),
        curve.move(pts.at(0)),
        ..pts.slice(1).map(p => curve.line(p))
      )
    )
  }
 })
]

= Output: txtual data representing scattering length

#set text(size: 10pt)
#for i in range(360) {
  [#i #fixed(cabs(psc.at(i)),6)]
  linebreak()
}

