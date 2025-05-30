import graph;

settings.outformat = "pdf";

size(500, 400, IgnoreAspect);  // Match other figures
scale(Log, Log);

// Load GeoPDEs data
string path1 = "geopdes_time.csv";
file f1 = input(path1);

real font = 17pt;

real[] dofs1, mean1;
bool first = true;
while (!eof(f1)) {
  string line = f1;
  if (first) { first = false; continue; }
  if (line == "") continue;
  string[] fields = split(line, ",");
  if (fields.length < 2) continue;
  dofs1.push((real)fields[0]);
  mean1.push((real)fields[1]);
}

// Load C++ solver data
string path2 = "cpp_time.csv";
file f2 = input(path2);

real[] dofs2, mean2;
first = true;
while (!eof(f2)) {
  string line = f2;
  if (first) { first = false; continue; }
  if (line == "") continue;
  string[] fields = split(line, ",");
  if (fields.length < 2) continue;
  dofs2.push((real)fields[0]);
  mean2.push((real)fields[1]);
}

// Compute log-space axis limits
real x_min = min(min(dofs1), min(dofs2));
real x_max = max(max(dofs1), max(dofs2));
real y_min = min(min(mean1), min(mean2));
real y_max = max(max(mean1), max(mean2));

xlimits(x_min / 1.5, x_max * 1.5);
ylimits(y_min / 3, y_max * 3);

// Plot GeoPDEs times
marker markG = marker(scale(1.4mm)*unitcircle, orange, Fill);
pen penG = orange + 2bp;
Label labG = Label("GeoPDEs", fontsize(font));
draw(graph(dofs1, mean1), penG, labG, markG);

// Plot C++ times
marker markC = marker(scale(1.4mm)*unitcircle, rgb(0, 0.5, 0), Fill);
pen penC = rgb(0, 0.5, 0) + 2bp;
Label labC = Label("C++ Solver", fontsize(font));
draw(graph(dofs2, mean2), penC, labC, markC);

// Axes and ticks
pen thin = invisible;
pen thick = gray + linetype("0 2") + linewidth(0.9);
pen minorTickPen = gray + 0.4bp;



Label xlabel = shift(0, -1.2)*Label("Degrees of Freedom", fontsize(font));
Label ylabel = shift(3mm*W)*rotate(90)*Label("Time (s)", fontsize(font));

xaxis(xlabel, BottomTop,
      LeftTicks(Label(fontsize(font)), begin=true, end=true, extend=true, ptick=thin, pTick=thick));
xaxis("", BottomTop,
      LeftTicks(format="%", ticklabel=null, ptick=minorTickPen, pTick=minorTickPen, extend=false));

yaxis(ylabel, LeftRight,
      RightTicks(Label(fontsize(font)), begin=true, end=true, extend=true, ptick=thin, pTick=thick));
yaxis("", LeftRight,
      RightTicks(format="%", ticklabel=null, ptick=minorTickPen, pTick=minorTickPen, extend=false));

// Legend and title
attach(legend(linelength=30bp, 1), point(SE), -20S + 20W, UnFill);
label(shift(2mm*N)*Label("\textbf{Diffusion problem}", fontsize(font)), point(N), N);