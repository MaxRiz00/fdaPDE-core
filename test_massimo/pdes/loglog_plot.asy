import graph;

size(450,350,IgnoreAspect);
scale(Log,Log);

// Load data
string path = "./results/sphere/L2_error.csv";
file f = input(path);

string header = f;

real[] x, y, z;
while (!eof(f)) {
  string line = f;
  if (line == "") continue;
  string[] fields = split(line, ",");
  if (fields.length < 3) continue;
  x.push((real)fields[0]);
  y.push((real)fields[1]);
  z.push((real)fields[2]); // H1 error
}

// Log-space padding
real x_factor = 1.5, y_factor = 3;
real x_min = min(x), x_max = max(x);
real y_min = min(min(y), min(z)), y_max = max(max(y), max(z));
xlimits(x_min / x_factor, x_max * x_factor);
ylimits(y_min / y_factor, y_max * y_factor);

// Plot L2 error
marker markL2 = marker(scale(1.2mm)*unitcircle, red, Fill);
pen dataPenL2 = red + 1.8bp;
Label err_label = Label("$\|u - u_h\|_{L^2(\Omega)}$", fontsize(15pt));
draw(graph(x, y), dataPenL2, err_label, markL2);

// Plot H1 error
marker markH1 = marker(scale(1.2mm)*unitcircle, blue, Fill);
pen dataPenH1 = blue + 1.8bp;
Label err_labelH1 = Label("$\|u - u_h\|_{H^1(\Omega)}$", fontsize(15pt));
draw(graph(x, z), dataPenH1, err_labelH1, markH1);

// Axes
pen thin = gray(0.7) + linetype("0 2") + linewidth(0.1);
pen thin2 = gray + linetype("0 2") + linewidth(0.9);
Label xlabel = shift(0,-1)*Label("$h$", fontsize(13pt));

xaxis(xlabel, BottomTop,
      LeftTicks(Label(fontsize(13pt)), begin=true, end=true, extend=true, ptick=thin, pTick=thin2));

yaxis(shift(2mm*W)*rotate(90)*"Error ", LeftRight,
      RightTicks(Label(fontsize(13pt)), begin=true, end=true, extend=true, ptick=thin, pTick=thin2));

// Reference slope lines
real x1 = 0.3, x2 = 0.09;
real[] refx = {x1, x2};

// h^3 line (L2 reference)
real y_ref = 2e-2;
real[] refy3 = {y_ref, y_ref * (refx[1]/refx[0])^3};
pen refPen3 = rgb(1, 0.6, 0.6) + linetype("4 2") + 1.8bp;
//Label h3_label = Label("$ h^3$", fontsize(15pt));
//draw(graph(refx, refy3), refPen3, h3_label );

// h^2 line (H1 reference)
real y_ref2 = 1e-1;
real[] refy2 = {y_ref2, y_ref2 * (refx[1]/refx[0])^2};
pen refPen2 = gray + linetype("4 2") + 1.8bp; //rgb(0.6, 0.6, 1) 
Label h2_label = Label("$ h^2$", fontsize(15pt));
draw(graph(refx, refy2), refPen2, h2_label );

// Attach legend
pen pleg = black;
attach(legend(linelength=10bp), point(SE), -15S + 13W, UnFill);

label(shift(2mm*N)*Label("\textbf{Sphere: } $p=2$", fontsize(16pt)), point(N), N);