import graph;

settings.outformat = "pdf";

size(600, 500, IgnoreAspect);
scale(Log, Log);

// Set font size and marker scale
real font = 21pt;
real markscale = 1.4mm;

// Data files and styles
string[] files = {
  "timings_nref_ring.csv",
  "timings_nref_curly.csv",
  "timings_nref_torus.csv",
  "timings_nref_twisted.csv"
};

string[] labels = {
  "Ring",
  "Curly plate",
  "Torus",
  "Twisted beam"
};

pen[] colors = {
  red,
  blue,
  darkgreen,
  orange
};

marker[] markers = {
  marker(scale(markscale)*unitcircle, red, Fill),
  marker(scale(markscale)*unitcircle, blue, Fill),
  marker(scale(markscale)*unitcircle, darkgreen, Fill),
  marker(scale(markscale)*unitcircle, orange, Fill)
};

// Read data
real[][] xs, ys;
for (int f = 0; f < files.length; ++f) {
  file file = input(files[f]);
  real[] x, y;
  while (!eof(file)) {
    string line = file;
    if (line == "" || line == "n_pts,time_sec") continue;
    string[] fields = split(line, ",");
    if (fields.length < 2) continue;
    x.push((real) fields[0]);
    y.push((real) fields[1]);
  }
  xs.push(x);
  ys.push(y);
}

// Set limits
real x_min = 1e10, x_max = -1e10, y_min = 1e10, y_max = -1e10;
for (int f = 0; f < xs.length; ++f) {
  x_min = min(x_min, min(xs[f]));
  x_max = max(x_max, max(xs[f]));
  y_min = min(y_min, min(ys[f]));
  y_max = max(y_max, max(ys[f]));
}
real x_factor = 1.5, y_factor = 1.5;
xlimits(x_min / x_factor, x_max * x_factor);
ylimits(y_min / y_factor, y_max * y_factor);

// Axes styling
pen thin = invisible;
pen thick = gray + linetype("0 2") + linewidth(0.9);
pen minorTickPen = gray + 0.4bp;

// x-axis
xaxis(shift(0, -1)*Label("$N_{cells}$", fontsize(font)), BottomTop,
      LeftTicks(Label(fontsize(font)), begin=true, end=true, extend=true, ptick=thin, pTick=thick));

xaxis("", BottomTop,
    LeftTicks(format = "%",
        begin=true, end=true, ticklabel=null,
        ptick=minorTickPen, pTick=minorTickPen, extend=false));

// y-axis
yaxis(shift(2.5mm*W)*rotate(90)*Label("Time (s)", fontsize(font)), LeftRight,
      RightTicks(Label(fontsize(font)), begin=true, end=true, extend=true, ptick=thin, pTick=thick));

yaxis("", LeftRight,
    RightTicks(format = "%",
        begin=true, end=true, ticklabel=null,
        ptick=minorTickPen, pTick=minorTickPen, extend=false));

// Plot all curves
for (int f = 0; f < files.length; ++f) {
  pen dataPen = colors[f] + 2bp;
  Label plotLabel = Label(labels[f], fontsize(font));
  draw(graph(xs[f], ys[f]), dataPen, plotLabel, markers[f]);
}

// Legend and title
attach(legend(linelength=30bp, 1), point(N), -35N + 40W, UnFill);
label(shift(2mm*N)*Label("\textbf{PI performance}", fontsize(font)), point(N), N);