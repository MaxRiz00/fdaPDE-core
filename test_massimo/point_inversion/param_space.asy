size(30cm);
import graph;
defaultpen(fontsize(10pt));

// === Define knot vectors ===
real[] xi = {0, 1};
real[] eta = {0, 0.111111, 0.222222, 0.333333, 0.444444, 0.555556, 0.666667, 0.777778, 0.888889, 1};

// === Draw knot lines (grid) ===
pen gridPen = gray + 0.6bp;

for (int i = 0; i < xi.length; ++i)
  draw((xi[i], eta[0]) -- (xi[i], eta[eta.length-1]), gridPen);

for (int i = 0; i < xi.length; ++i)
  label("$" + format(xi[i], "%.1f") + "$", (xi[i], -0.03), S, fontsize(4pt));

for (int j = 0; j < eta.length; ++j)
  draw((xi[0], eta[j]) -- (xi[xi.length-1], eta[j]), gridPen);

for (int j = 0; j < eta.length; ++j)
  label("$" + format(eta[j], "%.1f") + "$", (-0.02, eta[j]), W, fontsize(4pt));

// === Axes ===
draw((0,0)--(1.1,0), Arrow(6bp));
label("$\xi_1$", (1.1,0), E, fontsize(10pt));
draw((0,0)--(0,1.05), Arrow(6bp));
label("$\xi_2$", (0,1.05), N, fontsize(10pt));

// === Define descent path in parametric space ===
pair[] u = {
  (0.5, 0.888889),
  (0.3, 0.82141),
  (0.3, 0.847157),
  (0.3, 0.858246),
  (0.3, 0.86)
};

// === Plot arrows and dots ===
for (int i = 0; i < u.length - 1; ++i)
draw(u[i] -- u[i+1], black, Arrow(6bp));

for (int i = 0; i < u.length; ++i) {
  dot(u[i], red+3bp);
  if(i == 0) {
    label("$u_{0}$", u[i] + (0.01,0.01), NE, fontsize(5pt));
  }
  //label("$u_{" + string(i) + "}$", u[i] + (0.02,0.02), NE, fontsize(7pt));
}


// === Optional: mark final converged point ===
dot(u[u.length - 1], green+4bp);
label("$u_n \equiv u_{true}$", u[u.length - 1], N, fontsize(5pt));
