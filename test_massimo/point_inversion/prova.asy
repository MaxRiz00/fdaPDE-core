import graph;
size(8cm);
defaultpen(fontsize(10pt));

pair[] U = {
  (0.1, 0.3),
  (0.25, 0.4),
  (0.38, 0.47),
  (0.46, 0.51),
  (0.5, 0.5)
};

// Draw descent arrows
for (int i = 0; i < U.length - 1; ++i) {
  draw(U[i]--U[i+1], Arrow(6));
}

// Draw points
for (int i = 0; i < U.length; ++i) {
  dot(U[i], red);
  label("$u_{" + string(i) + "}$", U[i], NE);
}

/*
// Optional: Draw AABBs (example boxes around steps)
for (int i = 0; i < U.length; ++i) {
  pair p = U[i];
  draw(box(p - (0.05,0.05), p + (0.05,0.05)), gray+linewidth(0.5));
}
*/

// Optional: Plot a reference curve C(u)
real Cx(real u) { return u; }
real Cy(real u) { return 0.5*sin(2*pi*u) + 0.5; }

draw(graph(Cx, Cy, 0, 1, operator ..), dashed+blue);
label("Curve $C(u)$", (0.6, Cy(0.6)), SE, blue);